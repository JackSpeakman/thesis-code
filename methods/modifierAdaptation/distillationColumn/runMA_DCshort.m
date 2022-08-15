function [rk,uk,yk,objk,conk] = runMA_DCshort(varargin)
% run the distillation column using standard MA. 
% Both structrual and parametric uncertainty is present, with a reduced
% size of column and full efficient trays.

% Current results: converges slow ~11 iterates to G=0  (@K=0.6) or no
% convegrence (@K>=0.8)

%% 0. Deal with varargin
% default values
filter = 0.9;                        % filter gain
k_max = 21;                         % number of iterations
Ns = 14;                            % number of trays
Nr = 9;                             % feed tray
r_start = [70,85];                  % define starting point

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
    end
end

% print name
fprintf('\n#### Distilation Column (short) - %s - %4.2f - %2i/%2i ####\n\n','Standard MA',filter,Ns,Nr);

%% 1. Set-up variables
tic

% Define feed
Fin = [14/3600,0.32,71];
Fdiff = [32/30,30/32,1];

% Define parameters
N = Ns+Nr;
Nf = Ns;

theta_p = [8.5, 0.17, 0.155, 0.62,   0.35,   0.166, 0.5, 93900, 250, 190, 47.2];
theta_m = [8.5, 0.17, 0.155, 1, 1, 0.166, 0.51, 90900, 250, 190, 44.2, N, Nf];

% Set-up cost/constraint/controller functions
Profit = [200*32,400*60]*3.6;
cost = @(x,u)((u(2)*0.12));%-x(end)*Profit(2)-x(end-1)*Profit(1))+3.5);
cons = @(x,u,Nout)([x(Nout(1))-0.005,0.99-x(Nout(2))])*100;
Noutm = [3,N+3];
Noutp = [3,44];

contm = @(x)(x(2*N+2+[N-4,7]));
contp = @(x)(x(125+[28,14]));

tFun = @(t)(logspace(2,log10(t+100),10001)-100);
tspan = 60000000;

dr = diag([0.001,0.001]);

% Define plant
K = [-0.1,-0.0002,0.1,0.0002]; % Controller gain
M = diag([ones(1,84),zeros(1,83)]);
opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-8);
plant = @(r,x0)ode15s(@(t,x)plantPI(t,x,r,Fin.*Fdiff,K,theta_p),tFun(tspan),x0,opts);

% Define model
M = diag([ones(1,N+3),zeros(1,2*N+3)]);
opts = odeset('Mass',M,'RelTol',1e-10,'AbsTol',1e-10);
model = @(r,x0)ode15s(@(t,x)modelPI_short(t,x,r,Fin,K,theta_m),tFun(tspan),x0,opts);

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

fprintf('Set-up Complete       [%6.4fs]\n',toc)

%% 2. Find Initial Operating Point and Plant Optimum
tic
k = 1;

% Run model at some point
x0 = linspace(0.01,0.99,N+1);
V0 = linspace(0.0028,0.0032,N);
T0 = linspace(95,65,N+1);
Fout0 = [4/3600,10/3600]/69;

r0 = [71.7,84];
m0 = model(r0,[0,0,x0,V0,T0,Fout0]');
y0 = m0.y(:,end);

% Find model optimum
objF = @(r)(cOptFun(r,@(r)(model(r,y0)),@(r,y)(PIcontrol(r,K,y(1:2)',contm(y'))),@(u,y)(cost(y,u))));
conF = @(r)(cOptFun(r,@(r)(model(r,y0)),@(r,y)(PIcontrol(r,K,y(1:2)',contm(y'))),@(u,y)(cons(y,u,Noutm))));
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
fprintf('Beginning Standard MA run\n')
fprintf('%8s %10s %10s %10s %10s %10s %8s %10s [s]\n','k','T28','T14','Cost','Con 1','Con 2','Flag','Time');
fprintf('%8s %10.3f %10.3f %10.4f %10.4f %10.4f\n','init',rk(k,1),rk(k,2),objk(k),conk(k,1),conk(k,2))

for k = 2:k_max
    % calculate gradients
    m0 = model(rk(k-1,:),y0);
    y0 = m0.y(:,end);
    u0 = PIcontrol(rk(k-1,:),K,y0(1:2),contm(y0)).*[3600,1];
    obj0 = cost(y0,u0);
    con0 = cons(y0,u0,Noutm);
    
    dobjdr = zeros(2,1);
    dcondr = zeros(2,2);
    
    dobjpdr = zeros(2,1);
    dconpdr = zeros(2,2);
    
    for i = 1:2
        % shifted operating point
        r = rk(k-1,:) + dr(i,:);
        
        % model gradient
        m = model(r,y0);
        y = m.y(:,end);
        u = PIcontrol(r,K,y(1:2),contm(y)).*[3600,1];
        dobjdr(i) = (cost(y,u)-obj0)/dr(i,i);
        dcondr(i,:) = (cons(y,u,Noutm)-con0)/dr(i,i);
        
        % plant gradient
        p = plant(r,yp0);
        yp = p.y(:,end);
        up = PIcontrol(r,K,yp(1:2),contp(yp)).*[3600,1];
        dobjpdr(i) = (cost(yp,up)-objk(k-1,:))/dr(i,i);
        dconpdr(i,:) = (cons(yp,up,Noutp)-conk(k-1,:))/dr(i,i);
    end
    
    % calculate modified functions
    obj_mod0 = objk(k-1,:)-obj0;
    con_mod0 = conk(k-1,:)-con0;
    obj_mod1 = dobjpdr-dobjdr;
    con_mod1 = dconpdr-dcondr;
    
    objF = @(r)(cOptFun(r,@(r)(model(r,y0)),@(r,y)(PIcontrol(r,K,y(1:2)',contm(y')).*[3600,1]),@(u,y)(cost(y,u)))+...
        obj_mod0+(r-rk(k-1,:))*obj_mod1);    
    conF = @(r)(cOptFun(r,@(r)(model(r,y0)),@(r,y)(PIcontrol(r,K,y(1:2)',contm(y')).*[3600,1]),@(u,y)(cons(y,u,Noutm)))+...
        con_mod0+(r-rk(k-1,:))*con_mod1);

    % run optimization
    [rOpt,~,flag] = fmincon(@(r)objF(r),r0,[],[],[],[],[68,82],[74,86],@(r)(deal(conF(r),[])),opts);
    
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

fprintf('Finished Standard MA Run     [%6.4fs]\n',toc)

%% A. Embedded Functions
    function [c] = cOptFun(r,odeF,uF,cF)
        % Calculates the cost/constraints but does not re-run the
        % model/plant if the last run was for r (used in fmincon)
        
        %if any(r ~= r_last) 
            % run column
            r_last = r;
            ode = odeF(r_last);
            y_last = ode.y(:,end);
            u_last = uF(r_last, y_last);
        %end
        
        c = cF(u_last, y_last);
    end

    function resetHoldVar
        % results the embedded variables
        r_last = 0;
        y_last = 0;
        u_last = 0;
    end

end