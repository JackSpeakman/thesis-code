function [uk,yk,conk,objk] = runWCMA_WO(varargin)
% runs worst-case MA on the WO CSTR.
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
kmax = 31;          % Number of iterations
th = [0,0;70,160;-70,160;70,-160;-70,-160]*3;         % Default parameters
th_nom = [0,0];     % Default nominal parameters
Qk = 0;             % Use Qk or not (value of epsilon if yes)

conFun = @(u,y)WOconFun(u,y);   % constraint function

% replace certain values
n_in = floor(numel(varargin));
for i = 1:2:n_in
    if strcmp(varargin{i},'K')
        filterFun = varargin{i+1};
    elseif strcmp(varargin{i},'kmax')
        kmax = varargin{i+1};
    elseif strcmp(varargin{i},'th')
        th = varargin{i+1};
    elseif strcmp(varargin{i},'th_nom')
        th_nom = varargin{i+1};
    elseif strcmp(varargin{i},'conFun')
        conFun = varargin{i+1};
    elseif strcmp(varargin{i},'Qk')
        Qk = varargin{i+1};
    end
end

% replace K with a function
if isnumeric(filterFun)
    filterFun = @(~,~,~)(filterFun);
end

%% 1. Set-up parameters
% global variables
u_last = [0,0,0];
y_last = [0,0,0,0,0,0];

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
K = zeros(kmax,1);

%% 2. Find starting point
k = 1;
fminconopts = optimoptions('fmincon','Display','off');
fminimaxopts = optimoptions('fminimax','Display','off','MaxFunctionEvaluations',10000,'MaxIterations',1000);

uk(k,:) = fmincon(@(u)objFun(u,model(u,0)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,model(u,0)),[]),fminconopts);

yk(k,:) = plant(uk(k,:));
objk(k,:) = objFun(uk(k,:),yk(k,:));
conk(k,:) = conFun(uk(k,:),yk(k,:));

%% 3. Run MA
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
    
    xObj = @(u)(WCoptFun(u,model,objFun)+modO+(dobjpdu-dobjdu)'*(u-uk(k-1,:))');
    xCon = @(u)(WCoptFun(u,model,conFun)+modC+...
        permute(sum(repmat((u-uk(k-1,:))',1,n_th,n_c).*(dconpdu-dcondu),1),[2,3,1]));
    xQ = @(u)(quadFun((u-uk(k-1,:)),Qk,conk(k-1,:),dconpdu));
    
    % new optimum
    uOpt = fminimax(@(u)xObj(u),uGuess,[],[],[],[],umin,umax,...
        @(u)deal([xCon(u);xQ(u)],[]),fminimaxopts);
    
    % new operating point
    K(k) = filterFun(uk(k-1,:),uOpt,@(u)xCon(u)./norm(conk(k-1,:)));
    uk(k,:) = uk(k-1,:)*(1-K(k)) + uOpt*K(k);
    
    yk(k,:) = plant(uk(k,:));
    objk(k,:) = objFun(uk(k,:),yk(k,:));
    conk(k,:) = conFun(uk(k,:),yk(k,:));
    
end

%% 4. Embedded functions
    function out = WCoptFun(u,yFun,outFun)
        % calculates the optimization function for all th (efficient)
        if u == u_last
            y = y_last;
        else
            y = zeros(n_th,6);
            for jj = 1:n_th
                y(jj,:) = yFun(u,jj);
            end
            u_last = u;
            y_last = y;
        end
        out = outFun(u,y);
    end
end