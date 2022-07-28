function [uk,yk,conk,objk] = runPMAi_WO(varargin)
% runs probabilistic MA (individual) on the WO CSTR.
% ------------
% varargin          cell of inputs
%
% uk                inputs for iterations 1 to k
% yk                outputs for iterations 1 to k
% conk              constraints for iterations 1 to k
% objk              objective for iterations 1 to k
% ------------

%% 0. Deal with varargin
% default values
filter = 0.5;                       % Default filter function
kmax = 31;                          % Number of iterations
th_nom = [0,0];                     % Default nominal parameters
th = [0,0; 70,160; -70,160; 70,-160; -70,-160]*3;   % Default parameters
conFun = @(u,y)WOconFun(u,y);       % constraint function
u0 = [];                            % starting point
p = 0.9;                            % chance constraint [P(G_i(u)<=0) > p]

% replace certain values
n_in = floor(numel(varargin));
for i = 1:2:n_in
    switch varargin{i}
        case 'filter'
            filter = varargin{i+1};
        case 'kmax'
            kmax = varargin{i+1};
        case 'th'
            th = varargin{i+1};
        case 'th_nom'
            th_nom = varargin{i+1};
        case 'conFun'
            conFun = varargin{i+1};
        case 'startingPoint'
            u0 = varargin{i+1};
        case 'constraintChance'
            p = varargin{i+1};
    end
end

% print name
fprintf('\n#### Williams Otto - %s - %4.2f ####\n\n','Probabilistic MA (individual)',filter);

%% 1. Set-up parameters
tic

% global variables
u_last = [0,0,0];
y_last = [0,0,0,0,0,0];
thi_last = [0];

% WO functions
uGuess = [3.9,9.3,91];
yGuess = WOplantFun(uGuess);

n_th = size(th,1);
th = [th_nom;th];

model = @(u,i)WOmodelFun(u,yGuess,th(i+1,:));
plant = @(u,~)WOplantFun(u,yGuess);

objFun = @(u,y)WOobjFun(u,y);

% fix conFun size
size_c = size(conFun(uGuess,yGuess));
n_c = prod(size_c);
if size_c(1) ~= 1
    conFun = @(u,y)(conFun(u,y)');
end

% optimization limits
umin = [3,6,80];
umax = [4.5,11,105];
du = diag([0.0001,0.0001,0.0001]);

% proallocate outputs
uk = zeros(kmax,3);
yk = zeros(kmax,6);
objk = zeros(kmax,1);
conk = zeros(kmax,n_c);

%% 2. Find starting point
k = 1;
fminconopts = optimoptions('fmincon','Display','off');

% check if initial point is supplied
if isempty(u0)
    % if not, then solve model optimum
    uk(k,:) = fmincon(@(u)objFun(u,model(u,0)),uGuess,[],[],[],[],umin,umax,...
        @(u)deal(conFun(u,model(u,0)),[]),fminconopts);
else
    uk(k,:) = u0;
end

yk(k,:) = plant(uk(k,:));
objk(k,:) = objFun(uk(k,:),yk(k,:));
conk(k,:) = conFun(uk(k,:),yk(k,:));

%% 3. Run RTO
% print initial condition
t0 = toc;
fprintf('Beginning Worst Case MA run\n')
if n_c == 1
    fprintf('%8s %10s %10s %10s %10s %10s %8s %10s [s]\n','k','F_A','F_B','T_R','Cost','Con','Flag','Time');
    fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1))
elseif n_c == 2
    fprintf('%8s %10s %10s %10s %10s %10s %10s %8s %10s [s]\n','k','F_A','F_B','T_R','Cost','Con 1','Con 2','Flag','Time');
    fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f %10.4f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1),conk(k,2))
end

% run PMAi
for k = 2:kmax
    % gradients
    dobjpdu = zeros(3,1);
    dconpdu = zeros(3,1,n_c);
    dobjdu = zeros(3,n_th);
    dcondu = zeros(3,n_th,n_c);
    
    for i = 1:3
        u = uk(k-1,:) + du(i,:);
        
        ypi = plant(u);
        
        dobjpdu(i) = (objFun(u,ypi)-objFun(uk(k-1,:),yk(k-1,:)))/du(i,i);
        dconpdu(i,1,:) = (conFun(u,ypi)-conFun(uk(k-1,:),yk(k-1,:)))/du(i,i);
        for j = 1:n_th
            ym = model(uk(k-1,:),j);
            ymi = model(u,j);
            
            dobjdu(i,j) = (objFun(u,ymi)-objFun(uk(k-1,:),ym))/du(i,i);
            dcondu(i,j,:) = (conFun(u,ymi)-conFun(uk(k-1,:),ym))/du(i,i);
        end
    end
    
    % modified functions
    modO = zeros(n_th,1);
    modC = zeros(n_th,n_c);
    for j = 1:n_th
        modO(j,:) = objk(k-1,:) - objFun(uk(k-1,:),model(uk(k-1,:),j));
        modC(j,:) = conk(k-1,:) - conFun(uk(k-1,:),model(uk(k-1,:),j));
    end
    
    xObj = @(u)(WCoptFun(u,model,objFun,1:n_th)+modO+(dobjpdu-dobjdu)'*(u-uk(k-1,:))');
    xCon = @(u)(WCoptFun(u,model,conFun,1:n_th)+modC+...
        permute(sum(repmat((u-uk(k-1,:))',1,n_th,n_c).*(dconpdu-dcondu),1),[2,3,1]));
    
    % prob modified functions
    pObj = @(u)probCalc(xObj(u),0.5);
    pCon = @(u)probCalc(xCon(u),p);
    
    % new optimum
    [uOpt,~,flag] = fmincon(@(u)pObj(u),uGuess,[],[],[],[],umin,umax,...
        @(u)deal([pCon(u)],[]),fminconopts);
    
    % apply filter
    uk(k,:) = uk(k-1,:)*(1-filter) + uOpt*filter;
    
    % run plant
    yk(k,:) = plant(uk(k,:));
    objk(k,:) = objFun(uk(k,:),yk(k,:));
    conk(k,:) = conFun(uk(k,:),yk(k,:));
    
    % print result
    t1 = toc;
    if n_c == 1
        fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f %8i %10.3f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1),flag,t1-t0)
    elseif n_c == 2
        fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f %10.4f %8i %10.3f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1),conk(k,2),flag,t1-t0)
    end
    t0 = t1;
end

%% 4. Embedded functions
    function out = WCoptFun(u,yFun,outFun,thi)
        % calculates the optimization function for all th (efficient)
        nn_th = numel(thi);
        
        if all(u == u_last) && all(thi == thi_last)
            y = y_last;
        else
            y = zeros(nn_th,6);
            for jj = 1:nn_th
                y(jj,:) = yFun(u,thi(jj));
            end
            u_last = u;
            thi_last = thi;
            y_last = y;
        end
        out = outFun(u,y);
    end

    function pOut = probCalc(out,rho)
        % calculates the probabilistic value
        mu = mean(out);
        Sigma2 = var(out);
        zp = sqrt(2)*(erfinv(2*rho-1));
        pOut = mu+sqrt(Sigma2)*zp;
        
    end
end