function [uk,yk,conk,objk] = runWO2_MA(varargin)
% runs DHMA on the WO CSTR with 2 vars (fixed T) and parametric mismatch.
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
filterFun = 0.5;    % Default filter function
kmax = 21;          % Number of iterations
T = 90;             % fixed temperature


% replace certain values
n_in = floor(numel(varargin));
for i = 1:n_in
    if strcmp(varargin{i},'K')
        filterFun = varargin{i+1};
    elseif strcmp(varargin{i},'kmax')
        kmax = varargin{i+1};
    elseif strcmp(varargin{i},'T')
        T = varargin{i+1};
    end
end

% replace K with a function
if isnumeric(filterFun)
    filterFun = @(~,~,~)(filterFun);
    runK = 0;
else
    runK = 1;
end

%% 1. Set-up parameters
% global variables
u_last = [0,0];
y_last = [0,0,0,0,0,0];

% WO functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
uGuess = [3.88666741784971,9.36912123252326];

umin = [3,6];
umax = [4.5,11];

uNorm = @(u)((u-umin)./(umax-umin));
uRest = @(u)(umin+(umax-umin).*u);

th1 = [-200/3];
th2 = [200/3];
th3 = [2600/9];

th = [th1(:),th2(:),th3(:)];
n_th = size(th,1);

model = @(u,i)WOplantFun([uRest(u),T],yGuess,th(i+1,:));
plant = @(u)WOplantFun([uRest(u),T],yGuess);

objFun = @(u,y)WOobjFun([uRest(u),T],y);
conFun = @(u,y)WOconFun([uRest(u),T],y);

% fix conFun size
size_c = size(conFun(uGuess,yGuess));
n_c = prod(size_c);
if size_c(1) ~= 1
    conFun = @(u,y)(conFun(u,y)');
end

% optimization limits
umin = [3,6];
umax = [4.5,11];
du = diag([0.0001,0.0001]);

% proallocate outputs
uk = zeros(kmax,2);
yk = zeros(kmax,6);
objk = zeros(kmax,1);
conk = zeros(kmax,n_c);
K = zeros(kmax,1);

%% 2. Find starting point
k = 1;
K = 1;
fminopts = optimoptions('fmincon','Display','off','Algorithm','interior-point','MaxFunctionEvaluations',20000,'MaxIterations',10000);

uk(k,:) = fmincon(@(u)objFun(u,model(u,0)),uNorm(uGuess),[],[],[],[],[0,0],[1,1],...
    @(u)deal(conFun(u,model(u,0)),[]),fminopts);

uk(k,:) = [0.3,0.9];

yk(k,:) = plant(uk(k,:));
objk(k,:) = objFun(uk(k,:),yk(k,:));
conk(k,:) = conFun(uk(k,:),yk(k,:));

%% 3. Run MA
for k = 2:kmax
    % gradients    
    dobjpdu = zeros(2,1);
    dconpdu = zeros(2,1,n_c);
    dobjdu = zeros(2,1);
    dcondu = zeros(2,1,n_c);
    
    
    for i = 1:2
        u1 = uk(k-1,:) + du(i,:);
        ymi1 = model(u1,0);
        ypi1 = plant(u1);
        
        u2 = uk(k-1,:) - du(i,:);
        ymi2 = model(u2,0);
        ypi2 = plant(u2);
        
        dobjpdu(i) = (objFun(u1,ypi1)-objFun(u2,ypi2))/du(i,i)/2;
        dconpdu(i,1,:) = (conFun(u1,ypi1)-conFun(u2,ypi2))/du(i,i)/2;
        dobjdu(i) = (objFun(u1,ymi1)-objFun(u2,ymi2))/du(i,i)/2;
        dcondu(i,1,:) = (conFun(u1,ymi1)-conFun(u2,ymi2))/du(i,i)/2;
    end
    
    % modified functions
    u_last = 0;
    mod = conk(k-1,:) - conFun(uk(k-1,:),model(uk(k-1,:),0));
    
    xObj = @(u)(optFun(u,@(u)model(u,0),objFun)+(dobjpdu-dobjdu)'*(u-uk(k-1,:))');
    xCon = @(u)(optFun(u,@(u)model(u,0),conFun)+mod+(dconpdu-dcondu)'*(u-uk(k-1,:))');
    
    % new optimum
    uOpt = fmincon(@(u)xObj(u),uNorm(uGuess),[],[],[],[],[0,0],[1,1],...
        @(u)deal(xCon(u),[]),fminopts);
    
    % new operating point
    uk(k,:) = uk(k-1,:)*(1-K) + uOpt*K;
    
    yk(k,:) = plant(uk(k,:));
    objk(k,:) = objFun(uk(k,:),yk(k,:));
    conk(k,:) = conFun(uk(k,:),yk(k,:));
    
end

%% 4. Embedded functions
    function out = optFun(u,yFun,outFun)
        % calculates the optimization function (efficient)
        if u == u_last
            y = y_last;
        else
            y = yFun(u);
            u_last = u;
            y_last = y;
        end
        out = outFun(u,y);
    end
end