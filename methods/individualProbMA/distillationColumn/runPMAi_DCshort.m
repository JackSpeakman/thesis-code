function [rk,uk,yk,objk,conk] = runPMAi_DCshort(varargin)
% run the distillation column using probabilistic MA (individual). 
% Both structrual and parametric uncertainty is present, with a reduced
% size of column and full efficient trays.

% Current results: converges ~11 iterates to G=0 (K=0.9)

%% 0. Deal with varargin
% default values
filter = 0.9;                       % filter gain
k_max = 21;                         % number of iterations
Ns = 14;                            % number of trays
Nr = 7:11;                          % feed tray
r_start = [70,85];                  % define starting point
th_nomi = 3;                        % nominal rectifying tray
prob = 0.9;                         % chance constraint probability

for i = 1:2:nargin 
    switch varargin{i}
        case 'filter'
            filter = varargin{i+1};
        case 'k_max'
            k_max = varargin{i+1};
        case 'Ns'
            Ns = varargin{i+1};
        case 'Nr'
            Nr = varargin{i+1};
        case 'r_start'
            r_start = varargin{i+1};
        case 'th_nomi'
            th_nomi = varargin{i+1};
        case 'prob'
            prob = varargin{i+1};
    end
end

% print name
fprintf('\n#### Distilation Column (short) - %s - %4.2f - %2i/%2i ####\n\n','Probabilistic MA (individual)',filter,Ns,Nr);

%% 1. Set-up variables
tic

% Define feed
Fin = [14/3600,0.32,71];
Fdiff = [32/30,30/32,1];

% Define parameters (given to ode functions)
N = @(i)(Ns+Nr(i));     % total number of trays
Nf = Ns;                % feed tray
n_th = numel(Nr);       % number of models

theta_p = [8.5, 0.17, 0.155, 0.62,   0.35,   0.166, 0.5, 93900, 250, 190, 47.2];            % plant parameters
theta_m = @(i)([8.5, 0.17, 0.155, 1, 1, 0.166, 0.51, 90900, 250, 190, 44.2, N(i), Nf]);     % model parameters

% Set-up cost/constraint/controller functions
Profit = [200*32,400*60]*3.6;   % profit from distilate/bottoms

cost = @(x,u)((u(2)*0.12-x(end)*Profit(2)-x(end-1)*Profit(1))+3.7);
cons = @(x,u,Nout)([x(Nout(1))-0.005,0.986-x(Nout(2))])*100;

Noutp = [3,44];             % output trays (for plant)
Noutm = @(i)([3,N(i)+3]);   % output trays (for model)

contp = @(x)(x(125+[28,14]));               % function which gets key try temperatures from x_vector (plant)
contm = @(x,i)(x(2*N(i)+2+[N(i)-4,7]));     % function which gets key try temperatures from x_vector (model)

% time function (for ode solver)
tFun = @(t)(logspace(2,log10(t+100),10001)-100);
tspan = 600000000;

% step in r (for gradient calculation)
dr = diag([0.001,0.001]);

% Define plant
K = [-0.1,-0.0002,0.1,0.0002];              % Controller gain (plant and model)
M = diag([ones(1,84),zeros(1,83)]);         % ode mass matrix

opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-8);
plant = @(r,x0)ode15s(@(t,x)plantPI(t,x,r,Fin.*Fdiff,K,theta_p),tFun(tspan),x0,opts);

% Define model
M = @(i)(diag([ones(1,N(i)+3),zeros(1,2*N(i)+3)]));     % ode mass matrix

opts = @(i)(odeset('Mass',M(i),'RelTol',1e-10,'AbsTol',1e-10));
model = @(r,x0,i)ode15s(@(t,x)modelPI_short(t,x,r,Fin,K,theta_m(i)),tFun(tspan),x0,opts(i));

% Set-up variables for embedded functions used in the optimisation
r_last = 0;
y_last = 0;
u_last = 0;

% Allocate variables in RTO
objk = zeros(k_max,1);
conk = zeros(k_max,2);
yk = zeros(k_max,167);
uk = zeros(k_max,2);
rk = zeros(k_max,2);

% print 
fprintf('Set-up Complete       [%6.4fs]\n',toc)

%% 2. Find Initial Operating Point and Plant Optimum
tic
k = 1;

% Run all models at some point
y00 = cell(1,n_th);
for ith = 1:n_th
    % give some initial conditions to start ode solver from
    x0 = linspace(0.01,0.99,N(ith)+1);
    V0 = linspace(0.0028,0.0032,N(ith));
    T0 = linspace(95,65,N(ith)+1);
    Fout0 = [4/3600,10/3600]/69;
    
    % run model to get true column outputs
    r0 = [71.7,84];
    m0 = model(r0,[0,0,x0,V0,T0,Fout0]',ith);
    
    % save as initial column conditions
    y00{ith} = m0.y(:,end);
end

% Find model optimum
objF = @(r)(cOptFun(r,@(r)(model(r,y00{th_nomi},th_nomi)),@(r,y)(PIcontrol(r,K,y(1:2)',contm(y',th_nomi))),@(u,y)(cost(y,u))));
conF = @(r)(cOptFun(r,@(r)(model(r,y00{th_nomi},th_nomi)),@(r,y)(PIcontrol(r,K,y(1:2)',contm(y',th_nomi))),@(u,y)(cons(y,u,Noutm(i)))));
opts = optimoptions('fmincon','Display','off','Algorithm','sqp','StepTolerance',1e-9,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);

if r_start(1) == 0
    rk(k,:) = fmincon(@(r)objF(r),r0,[],[],[],[],[68,82],[74,88],@(r)(deal(conF(r)+[0.1,0],[])),opts);
else
    rk(k,:) = r_start;
end

% Run plant at some point
x0 = [linspace(0.001,0.01,5),linspace(0.012,0.988,32),linspace(0.99,0.999,5)];
n0 = linspace(0.0034,0.0044,40);
V0 = linspace(0.0028,0.0032,40);
T0 = linspace(95,65,41);
Fout0 = [4/3600,10/3600]/69;

r0 = [71.7,84];
p0 = plant(r0,[0,0,x0,n0,V0,T0,Fout0]');
yp0 = p0.y(:,end);

% Find plant optimum
resetHoldVar
rO = 0;

if rO == 1
    objF = @(r)(cOptFun(r,@(r)(plant(r,yp0)),@(r,y)(PIcontrol(r,K,y(1:2)',contp(y')).*[3600,1]),@(u,y)(cost(y,u))));
    conF = @(r)(cOptFun(r,@(r)(plant(r,yp0)),@(r,y)(PIcontrol(r,K,y(1:2)',contp(y')).*[3600,1]),@(u,y)(cons(y,u,Noutp))));
    
    rp_Opt = fmincon(@(r)objF(r),r0,[],[],[],[],[68,82],[74,86],@(r)(deal(conF(r),[])),opts);
else
    rp_Opt = [69.8863632317303,83.2430177445115];
end

% Run plant at starting point
resetHoldVar
p = plant(rk(k,:),yp0);
yk(k,:) = p.y(:,end);
uk(k,:) = PIcontrol(rk(k,:),K,yk(k,1:2)',contp(yk(k,:)')).*[3600,1];
objk(k,:) = cost(yk(k,:),uk(k,:));
conk(k,:) = cons(yk(k,:),uk(k,:),Noutp);

fprintf('Plant Optimum Found   [%6.4fs]\n',toc),tic

%% 3. Run Standard MA
t0 = toc;
fprintf('Beginning Worst-Case MA run\n')
fprintf('%8s %10s %10s %10s %10s %10s %8s %10s [s]\n','k','T28','T14','Cost','Con 1','Con 2','Flag','Time');
fprintf('%8s %10.3f %10.3f %10.4f %10.4f %10.4f\n','init',rk(k,1),rk(k,2),objk(k),conk(k,1),conk(k,2))

for k = 2:k_max
    % pre allocate
    dobjdr = zeros(2,1,n_th);
    dcondr = zeros(2,2,n_th);
    dobjpdr = zeros(2,1);
    dconpdr = zeros(2,2);
    obj0 = zeros(n_th,1);
    con0 = zeros(n_th,2);
    
    % plant gradient
    for i = 1:2
        % positive shifted operating point
        r1 = rk(k-1,:) + dr(i,:);
        p1 = plant(r1,yp0);
        yp1 = p1.y(:,end);
        up1 = PIcontrol(r1,K,yp1(1:2),contp(yp1)).*[3600,1];
        
        % negative shifted operating point
        r2 = rk(k-1,:) - dr(i,:);
        p2 = plant(r2,yp0);
        yp2 = p2.y(:,end);
        up2 = PIcontrol(r2,K,yp2(1:2),contp(yp2)).*[3600,1];
        
        % calc gradient
        dobjpdr(i) = (cost(yp1,up1)-cost(yp2,up2))/(2*dr(i,i));
        dconpdr(i,:) = (cons(yp1,up1,Noutp)-cons(yp2,up2,Noutp))/(2*dr(i,i));
    end
            
    for ith = 1:n_th
        % calculate gradients
        m0 = model(rk(k-1,:),y00{ith},ith);
        y0 = m0.y(:,end);
        u0 = PIcontrol(rk(k-1,:),K,y0(1:2),contm(y0,ith)).*[3600,1];
        obj0(ith) = cost(y0,u0);
        con0(ith,:) = cons(y0,u0,Noutm(ith));
        
        for i = 1:2
            % positive shifted operating point
            r1 = rk(k-1,:) + dr(i,:);
            m1 = model(r1,y0,ith);
            y1 = m1.y(:,end);
            u1 = PIcontrol(r1,K,y1(1:2),contm(y1,ith)).*[3600,1];
            
            % negative shifted operating point
            r2 = rk(k-1,:) - dr(i,:);
            m2 = model(r2,y0,ith);
            y2 = m2.y(:,end);
            u2 = PIcontrol(r2,K,y2(1:2),contm(y2,ith)).*[3600,1];
            
            % model gradient
            dobjdr(i,ith) = (cost(y1,u1)-cost(y2,u2))/(2*dr(i,i));
            dcondr(i,:,ith) = (cons(y1,u1,Noutm(ith))-cons(y2,u2,Noutm(ith)))/(2*dr(i,i));
            
        end
    end
    % calculate modified functions
    obj_mod0 = objk(k-1,:)-obj0;
    con_mod0 = conk(k-1,:)-con0;
    obj_mod1 = dobjpdr-dobjdr;
    con_mod1 = dconpdr-dcondr;
    
    objF = @(r)(cOptFunAll(r,@(r,i)(model(r,y00{i},i)),@(r,y,i)(PIcontrol(r,K,y(1:2)',contm(y',i)).*[3600,1]),@(u,y,i)(cost(y,u)),1)+...
        obj_mod0+sum(repmat((r-rk(k-1,:)),5,1).*permute(obj_mod1,[3,1,2]),2));  
    
    conF = @(r)(cOptFunAll(r,@(r,i)(model(r,y00{i},i)),@(r,y,i)(PIcontrol(r,K,y(1:2)',contm(y',i)).*[3600,1]),@(u,y,i)(cons(y,u,Noutm(i))),2)+...
        con_mod0+permute(sum(repmat((r-rk(k-1,:)),5,1,2).*permute(con_mod1,[3,1,2]),2),[1,3,2]));

    % prob modified functions
    pObj = @(u)probCalc(objF(u),0.5);
    pCon = @(u)probCalc(conF(u),prob);
    
    % run optimization
    [rOpt,~,flag] = fmincon(@(r)pObj(r),r0,[],[],[],[],[68,82],[74,86],@(r)(deal(pCon(r),[])),opts);
    
    % run plant at new operating point
    rk(k,:) = rOpt*filter + rk(k-1,:)*(1-filter);
    
    p = plant(rk(k,:),yp0);
    yk(k,:) = p.y(:,end);
    uk(k,:) = PIcontrol(rk(k,:),K,yk(k,1:2)',contp(yk(k,:)')).*[3600,1];
    objk(k,:) = cost(yk(k,:),uk(k,:));
    conk(k,:) = cons(yk(k,:),uk(k,:),Noutp);
    
    % print results
    t1 = toc;
    fprintf('%8i %10.3f %10.3f %10.4f %10.4f %10.4f %8i %10.3f\n',k,rk(k,1),rk(k,2),objk(k),conk(k,1),conk(k,2),flag,t1-t0)
    t0 = t1;
end

fprintf('Finished PMAi Run     [%6.4fs]\n',toc)

%% A. Embedded Functions
    function [c] = cOptFun(r,odeF,uF,cF)
        % Calculates the cost/constraints but does not re-run the
        % model/plant if the last run was for r (used in fmincon)
        
        if any(r ~= r_last) 
            % run column
            r_last = r;
            ode = odeF(r_last);
            y_last = ode.y(:,end);
            u_last = uF(r_last,y_last);
        end
        
        c = cF(u_last,y_last);
    end

    function [c] = cOptFunAll(r,odeF,uF,cF,n_c)
        % Calculates the cost/constraints but does not re-run the
        % model/plant if the last run was for r (used in fmincon)
        
       if any(r ~= r_last) 
            % run column
            r_last = r;
            y_last = cell(1,n_th);
            u_last = cell(1,n_th);
            
            for iith = 1:n_th
                ode = odeF(r_last,iith);
                y_last{iith} = ode.y(:,end);
                u_last{iith} = uF(r_last,y_last{iith},iith);
            end
       end
        
        c = zeros(n_th,n_c);
        for iith = 1:n_th
            c(iith,:) = cF(u_last{iith},y_last{iith},iith);
        end
    end

    function pOut = probCalc(out,rho)
        % calculates the probabilistic value
        mu = mean(out);
        Sigma2 = var(out);
        zp = sqrt(2)*(erfinv(2*rho-1));
        pOut = mu+sqrt(Sigma2)*zp;
    end

    function resetHoldVar
        % results the embedded variables
        r_last = 0;
        y_last = 0;
        u_last = 0;
    end

end