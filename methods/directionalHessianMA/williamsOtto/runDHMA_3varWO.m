function [uk,yk,conk,objk] = runDHMA_3varWO(varargin)
% runs DHMA on the WO CSTR with 3 vars (fixed T) and parametric mismatch.
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
K = 0.5;    % Default filter function
kmax = 21;          % Number of iterations
fileName = 'Sbar\SbarData\Sbar03.mat'; % file name


% replace certain values
n_in = floor(numel(varargin));
for i = 1:n_in
    if strcmp(varargin{i},'kmax')
        kmax = varargin{i+1};
    elseif strcmp(varargin{i},'fileName')
        fileName = varargin{i+1};
    elseif strcmp(varargin{i},'K')
        K = varargin{i+1};
    end
end

% load data
data = load(fileName);

%% 1. Set-up parameters
% WO functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
uGuess = [3.88666741784971,9.36912123252326,91];

umin = [3,6,80];
umax = [4.5,11,105];

uNorm = @(u)((u-umin)./(umax-umin));
uRest = @(u)(umin+(umax-umin).*u);

model = @(u)WOmodelFun(uRest(u),yGuess);
plant = @(u)WOplantFun(uRest(u),yGuess);

objFun = @(u,y)WOobjFun(uRest(u),y);
conFun = @(u,y)WOconFun(uRest(u),y);

% optimization limits
du = diag([0.001,0.001,0.001]);

% proallocate outputs
uk = zeros(kmax,3);
yk = zeros(kmax,6);
objk = zeros(kmax,1);
conk = zeros(kmax,1);

%% 2. Find starting point
k = 1;
fminopts = optimoptions('fmincon','Display','off','Algorithm','interior-point','MaxFunctionEvaluations',20000,'MaxIterations',10000);

uk(k,:) = fmincon(@(u)objFun(u,model(u)),uNorm(uGuess),[],[],[],[],[0,0,0],[1,1,1],...
    @(u)deal(conFun(u,model(u)),[]),fminopts);

yk(k,:) = plant(uk(k,:));
objk(k,:) = objFun(uk(k,:),yk(k,:));
conk(k,:) = conFun(uk(k,:),yk(k,:));

%% 3. Run DHMA
for k = 2:kmax
    % gradients    
    dobjpdu = zeros(3,1);
    dconpdu = zeros(3,1);
    
    for i = 1:3
        u1 = uk(k-1,:) + du(i,:);
        ypi1 = plant(u1);
        
        u2 = uk(k-1,:) - du(i,:);
        ypi2 = plant(u2);
        
        dobjpdu(i) = (objFun(u1,ypi1)-objFun(u2,ypi2))/du(i,i)/2;
        dconpdu(i) = (conFun(u1,ypi1)-conFun(u2,ypi2))/du(i,i)/2;
    end
    
    % modified functions
    xObj = @(u)(objk(k-1,:)+(dobjpdu)'*(u-uk(k-1,:))'+calcDelta3D(uk(k-1,:),u,data.vv1,data.vv2,data.vv3,data.az,data.el,data.Tri,data.Sbar.obj));
    xCon = @(u)(conk(k-1,:)+(dconpdu)'*(u-uk(k-1,:))'+calcDelta3D(uk(k-1,:),u,data.vv1,data.vv2,data.vv3,data.az,data.el,data.Tri,data.Sbar.con));
    
    % new optimum
    uOpt = fmincon(@(u)xObj(u),uNorm(uGuess),[],[],[],[],[0,0,0],[1,1,1],...
        @(u)deal(xCon(u),[]),fminopts);
    
    % new operating point
    uk(k,:) = uk(k-1,:)*(1-K) + uOpt*K;
    
    yk(k,:) = plant(uk(k,:));
    objk(k,:) = objFun(uk(k,:),yk(k,:));
    conk(k,:) = conFun(uk(k,:),yk(k,:));
    
end

%% 4. Embedded functions

end