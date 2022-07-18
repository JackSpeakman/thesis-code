function [uk,yk,conk,objk] = runMA_WO(varargin)
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
% guess variables
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
uGuess = [3.88, 9.36, 91.23];

% default values
filter = 0.5;                       % Default filter function
kmax = 21;                          % Number of iterations
conFun = @(u,y)WOconFun(u,y);       % constraint function
u0 = [];                            % starting point
model = @(u)WOmodelFun(u,yGuess);   % model function
plant = @(u)WOplantFun(u,yGuess);   % plant function
n_u = 3;                            % number of inputs

% replace certain values
n_in = floor(numel(varargin));
for i = 1:2:n_in
    switch varargin{i}
        case 'filter'
            filter = varargin{i+1};
        case 'kmax'
            kmax = varargin{i+1};
        case 'conFun'
            conFun = varargin{i+1};
        case 'startingPoint'
            u0 = varargin{i+1};
        case 'modelFun'
            model = varargin{i+1};
        case 'plantFun'
            plant = varargin{i+1};
        case 'num_inputs'
            n_u = varargin{i+1};
    end
end

% print name
fprintf('\n#### Williams Otto - %ivar - %s - %4.2f ####\n\n',n_u,'Standard MA',filter);

%% 1. Set-up parameters
tic

% global variables
u_last = zeros(1,n_u);
y_last = [0,0,0,0,0,0];

% objective functions
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
du = diag(ones(1,n_u)*0.0001);

% proallocate outputs
uk = zeros(kmax,n_u);
yk = zeros(kmax,6);
objk = zeros(kmax,1);
conk = zeros(kmax,n_c);
K = zeros(kmax,1);


%% 2. Find starting point
k = 1;
fminopts = optimoptions('fmincon','Display','off','Algorithm','interior-point','MaxFunctionEvaluations',20000,'MaxIterations',10000);

% check if initial point is supplied
if isempty(u0)
    % if not, then solve model optimum
    uk(k,:) = fmincon(@(u)objFun(u,model(u)),uGuess(1:n_u),[],[],[],[],umin(1:n_u),umax(1:n_u),...
        @(u)deal(conFun(u,model(u)),[]),fminopts);
else
    uk(k,:) = u0;
end

% run plant
yk(k,:) = plant(uk(k,:));
objk(k,:) = objFun(uk(k,:),yk(k,:));
conk(k,:) = conFun(uk(k,:),yk(k,:));

fprintf('Set-up Complete       [%6.4fs]\n',toc)

%% 3. Run MA
% print initial condition
t0 = toc;
fprintf('Beginning Standard MA run\n')
if n_u == 2 && n_c == 1
    fprintf('%8s %10s %10s %10s %10s %8s %10s [s]\n','k','F_A','F_B','Cost','Con','Flag','Time');
    fprintf('%8i %10.3f %10.3f %10.4f %10.4f\n',k,uk(k,1),uk(k,2),objk(k),conk(k,1))
elseif n_u == 3 && n_c == 1
    fprintf('%8s %10s %10s %10s %10s %10s %8s %10s [s]\n','k','F_A','F_B','T_R','Cost','Con','Flag','Time');
    fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1))
elseif n_u == 2 && n_c == 2
    fprintf('%8s %10s %10s %10s %10s %10s %8s %10s [s]\n','k','F_A','F_B','Cost','Con 1','Con 2','Flag','Time');
    fprintf('%8i %10.3f %10.3f %10.4f %10.4f %10.4f\n',k,uk(k,1),uk(k,2),objk(k),conk(k,1),conk(k,2))
elseif n_u == 3 && n_c == 2
    fprintf('%8s %10s %10s %10s %10s %10s %10s %8s %10s [s]\n','k','F_A','F_B','T_R','Cost','Con 1','Con 2','Flag','Time');
    fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f %10.4f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1),conk(k,2))
end

% run MA
for k = 2:kmax
    % pre-allocate gradients    
    dobjpdu = zeros(n_u,1);
    dconpdu = zeros(n_u,n_c);
    dobjdu = zeros(n_u,1);
    dcondu = zeros(n_u,n_c);
    
    % model at initial point
    for i = 1:n_u
        % positive shift
        u_pos = uk(k-1,:) + du(i,:);
        ymi_pos = model(u_pos);
        ypi_pos = plant(u_pos);
        
        % negative shift
        u_neg = uk(k-1,:) - du(i,:);
        ymi_neg = model(u_neg);
        ypi_neg = plant(u_neg);
        
        % get gradients
        dobjpdu(i) =   (objFun(u_pos,ypi_pos)-objFun(u_neg,ypi_neg))/(2*du(i,i));
        dconpdu(i,:) = (conFun(u_pos,ypi_pos)-conFun(u_neg,ypi_neg))/(2*du(i,i));
        
        dobjdu(i) =   (objFun(u_pos,ymi_pos)-objFun(u_neg,ymi_neg))/(2*du(i,i));
        dcondu(i,:) = (conFun(u_pos,ymi_pos)-conFun(u_neg,ymi_neg))/(2*du(i,i));
    end
    
    % modified functions
    u_last = 0;
    mod = conk(k-1,:) - conFun(uk(k-1,:),model(uk(k-1,:)));
    
    xObj = @(u)(optFun(u,@(u)model(u),objFun)+(u-uk(k-1,:))*(dobjpdu-dobjdu));
    xCon = @(u)(optFun(u,@(u)model(u),conFun)+(u-uk(k-1,:))*(dconpdu-dcondu)+mod);
    
    % new optimum
    [uOpt,~,flag] = fmincon(@(u)xObj(u),uGuess(1:n_u),[],[],[],[],umin(1:n_u),umax(1:n_u),...
        @(u)deal(xCon(u),[]),fminopts);
    
    % apply filter
    uk(k,:) = uk(k-1,:)*(1-filter) + uOpt*filter;
    
    % run plant
    yk(k,:) = plant(uk(k,:));
    objk(k,:) = objFun(uk(k,:),yk(k,:));
    conk(k,:) = conFun(uk(k,:),yk(k,:));
    
    % print result
    t1 = toc;
    if n_u == 2 && n_c == 1
        fprintf('%8i %10.3f %10.3f %10.4f %10.4f %8i %10.3f\n',k,uk(k,1),uk(k,2),objk(k),conk(k,1),flag,t1-t0)
    elseif n_u == 3 && n_c == 1
        fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f %8i %10.3f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1),flag,t1-t0)
    elseif n_u == 2 && n_c == 2
        fprintf('%8i %10.3f %10.3f %10.4f %10.4f %10.4f %8i %10.3f\n',k,uk(k,1),uk(k,2),objk(k),conk(k,1),conk(k,2),flag,t1-t0)
    elseif n_u == 3 && n_c == 2
        fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f %10.4f %8i %10.3f\n',k,uk(k,1),uk(k,2),uk(k,3),objk(k),conk(k,1),conk(k,2),flag,t1-t0)
    end
    t0 = t1;
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