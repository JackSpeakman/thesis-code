function [uk,yk,conk,objk] = runWO_MA(varargin)
% runs standard MA on the WO CSTR.
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
th = [0,0;70,160;-70,160;70,-160;-70,-160]*3;   % Default parameters (used by filter)
th_nom = [0,0];     % Default nominal parameters
Qk = 0;             % Use Qk or not (value of epsilon if yes)

conFun = @(u,y)WOconFun(u,y);   % constraint function

% replace certain values
n_in = floor(numel(varargin));
for i = 1:n_in
    if strcmp(varargin{i},'K')
        filterFun = varargin{i+1};
    elseif strcmp(varargin{i},'kmax')
        kmax = varargin{i+1};
    elseif strcmp(varargin{i},'th')
        th = varargin{i+1};
    elseif strcmp(varargin{i},'th_nom')
        th = varargin{i+1};
    elseif strcmp(varargin{i},'conFun')
        conFun = varargin{i+1};
    elseif strcmp(varargin{i},'Qk')
        Qk = varargin{i+1};
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
u_last = [0,0,0];
y_last = [0,0,0,0,0,0];

% WO functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
uGuess = [3.88666741784971,9.36912123252326,91.2326607595759];

n_th = size(th,1);
th = [th_nom;th];

model = @(u,i)WOmodelFun(u,yGuess,th(i+1,:));
plant = @(u)WOplantFun(u,yGuess);

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
du = diag([0.001,0.001,0.001]);

% proallocate outputs
uk = zeros(kmax,3);
yk = zeros(kmax,6);
objk = zeros(kmax,1);
conk = zeros(kmax,n_c);
K = zeros(kmax,1);

%% 2. Find starting point
k = 1;
fminopts = optimoptions('fmincon','Display','final','Algorithm','interior-point','MaxFunctionEvaluations',20000,'MaxIterations',10000);

uk(k,:) = fmincon(@(u)objFun(u,model(u,0)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,model(u,0)),[]),fminopts);

yk(k,:) = plant(uk(k,:));
objk(k,:) = objFun(uk(k,:),yk(k,:));
conk(k,:) = conFun(uk(k,:),yk(k,:));

%% 3. Run MA
for k = 2:kmax
    % gradients    
    dobjpdu = zeros(3,1);
    dconpdu = zeros(3,1,n_c);
    dobjdu = zeros(3,1);
    dcondu = zeros(3,1,n_c);
    
    ym = model(uk(k-1,:),0);
    
    for i = 1:3
        u = uk(k-1,:) + du(i,:);
        
        ymi = model(u,0);
        ypi = plant(u);
        
        dobjpdu(i) = (objFun(u,ypi)-objFun(uk(k-1,:),yk(k-1,:)))/du(i,i);
        dconpdu(i,1,:) = (conFun(u,ypi)-conFun(uk(k-1,:),yk(k-1,:)))/du(i,i);
        dobjdu(i) = (objFun(u,ymi)-objFun(uk(k-1,:),ym))/du(i,i);
        dcondu(i,1,:) = (conFun(u,ymi)-conFun(uk(k-1,:),ym))/du(i,i);
    end
    
    % modified functions
    u_last = 0;
    mod = conk(k-1,:) - conFun(uk(k-1,:),model(uk(k-1,:),0));
    
    xObj = @(u)(optFun(u,@(u)model(u,0),objFun)+(dobjpdu-dobjdu)'*(u-uk(k-1,:))');
    xCon = @(u)(optFun(u,@(u)model(u,0),conFun)+mod+...
        permute(sum(repmat((u-uk(k-1,:))',1,n_th,n_c).*(dconpdu-dcondu),1),[2,3,1]));
    xQ = @(u)(quadFun((u-uk(k-1,:)),Qk,conk(k-1,:),dconpdu));
    
    % new optimum
    uOpt = fmincon(@(u)xObj(u),uGuess,[],[],[],[],umin,umax,...
        @(u)deal([xCon(u);xQ(u)],[]),fminopts);
    
    % new operating point
    if runK
        % calculate modified function for the filter
        modK = zeros(n_th,n_c);
        dconduK = zeros(3,n_th,n_c);
        dconpduK = zeros(3,1,n_c);
        for i = 1:3
            u = uk(k-1,:) + du(i,:);
            ypi = plant(u);
            
            dconpduK(i,1,:) = (conFun(u,ypi)-conFun(uk(k-1,:),yk(k-1,:)))/du(i,i);
            for j = 1:n_th
                ym = model(uk(k-1,:),j);
                ymi = model(u,j);
                
                dconduK(i,j,:) = (conFun(u,ymi)-conFun(uk(k-1,:),ym))/du(i,i);
            end
        end
        
        for j = 1:n_th
            modK(j,:) = conk(k-1,:) - conFun(uk(k-1,:),model(uk(k-1,:),j));
        end
        u_last = 0;
        xConK = @(u)(WCoptFun(u,model,conFun)+modK+...
            permute(sum(repmat((u-uk(k-1,:))',1,n_th,n_c).*(dconpduK-dconduK),1),[2,3,1]));
    end
        
    K(k) = filterFun(uk(k-1,:),uOpt,@(u)xConK(u)./norm(conk(k-1,:)));
    uk(k,:) = uk(k-1,:)*(1-K(k)) + uOpt*K(k);
    
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