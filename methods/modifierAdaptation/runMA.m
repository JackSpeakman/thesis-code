function [uk,yk,conk,objk] = runMA(varargin)
% The master code for running the modifier adaptation method. The standard
% MA approach is described in the I&EC article titled "Modifier-Adaptation 
% Methodology for Real-Time Optimization" by A. Marchetti, B. Chachuat, and 
% D. Bonvin. 
% 
% ------ INPUT VARIABLES ------
% varargin                                  cell of inputs
%   'filter'            1-by-1              Input filter
%   'kmax'              1-by-1              Number of RTO iterations
%   'stratingPoint'     1-by-n_u            RTO starting point
%   'uGuess'            1-by-n_u            Input guess (for initial runs)
%   'conFun'            @(u,y)              Constraint function
%   'objFun'            @(u,y)              Objective function
%   'modelFun'          @(u,th)             Model function
% 
% ------ OUTPUT VARIABLES ------
% uk        kmax-by-n_u         Inputs for iterations 1 to kmax
% yk        kmax-by-n_y         Outputs for iterations 1 to kmax
% conk      kmax-by-n_c         Constraints for iterations 1 to kmax
% objk      kmax-by-1           Objective for iterations 1 to kmax
% 
% ------ EXAMPLES ------
% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runMA
%       This runs the default argument which is MA applied to WO
%
% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runRobustMA('filter',0.7,'StartingPoint',[3.2,6.5,90])
%       This runs MA applied to WO with an input filter gain of 0.7, 
%       and the initial point of [3.2,6.5,90]
% 
% 
% See 'distillationColumn/makePlots___.m' for distillation colunm examples
% See 'semiBatchReactor/makePlots___.m' for semi-batch reactor examples

%% 0. Deal with varargin
% set-up default values [for MA on WO]
% general RTO parameters
filter = 0.5;                       % input filter
kmax = 11;                          % number of iterations
u0 = [];                            % starting point

% model/plant parameters
conFun = @(u,y)WOconFun(u,y);       % constraint function
objFun = @(u,y)WOobjFun(u,y);       % objective function
model = @(u,yG)WOmodelFun(u,yG);    % model function
plant = @(u)(WOplantFun(u));        % plant function
umin = [3,6,80];                    % optimization limit [min]
umax = [4.5,11,105];                % optimization limit [max]

% replace certain values
n_in = floor(numel(varargin));
for i = 1:2:n_in
    switch varargin{i}
        case 'filter'
            filter = varargin{i+1};
        case 'kmax'
            kmax = varargin{i+1};
        case 'startingPoint'
            u0 = varargin{i+1};
        case 'uGuess'
            uGuess = varargin{i+1};
        case 'th'
            th = varargin{i+1};
        case 'th_nom'
            th_nom = varargin{i+1};
        case 'conFun'
            conFun = varargin{i+1};
        case 'objFun'
            conFun = varargin{i+1};
        case 'modelFun'
            model_th = varargin{i+1};
        case 'method'
            method = varargin{i+1};
        case 'umin'
            umin = varargin{i+1};
        case 'umax'
            umax = varargin{i+1};
        case 'probChance'
            rho = varargin{i+1};
    end
end

% print name
fprintf('\n#### Standard MA - %4.2f ####\n\n',filter);

%% 1. Set-up parameters
tic

% uGuess
if ~exist('uGuess','var')
    uGuess = (umin+umax)/2;
end
n_u = numel(uGuess);

% get yGuess
yGuess = model(uGuess,[]);
model_u = @(u)model(u,yGuess);

% fix conFun size
size_c = size(conFun(uGuess,yGuess));
n_c = prod(size_c);
if size_c(1) ~= 1
    conFun = @(u,y)(conFun(u,y)');
end

% get model outputs size
n_y = numel(yGuess);

% gradient shift
du = eye(n_u)*0.0001;

% proallocate outputs
uk = zeros(kmax,n_u);
yk = zeros(kmax,n_y);
objk = zeros(kmax,1);
conk = zeros(kmax,n_c);

% global variables
u_last = zeros(1,n_u);
y_last = zeros(1,n_y);

%% 2. Find starting point
k = 1;
fminopts = optimoptions('fmincon','Display','off');

% check if initial point is supplied
if isempty(u0)
    % if not, then solve model optimum
    uk(k,:) = fmincon(@(u)objFun(u,model_u(u)),uGuess,[],[],[],[],umin,umax,...
        @(u)deal(conFun(u,model_u(u)),[]),fminopts);
else
    uk(k,:) = u0;
end

yk(k,:) = plant(uk(k,:));
objk(k,:) = objFun(uk(k,:),yk(k,:));
conk(k,:) = conFun(uk(k,:),yk(k,:));

fprintf('Set-up Complete       [%6.4fs]\n',toc)

%% 3. Run MA
% print initial condition
t0 = toc;
fprintf('Beginning Standard MA run\n')
fprintf(['%8s ' repmat('%10s ',1,n_u) '%10s ' repmat('%10s ',1,n_c) '%8s %10s [s]\n'],...
    'k',string([repmat('Input ',n_u,1),num2str((1:n_u)')]),'Cost',string([repmat('Con ',n_c,1),num2str((1:n_c)')]),'Flag','Time');

fprintf('%8s %10.3f %10.3f %10.3f %10.4f %10.4f\n','init',uk(k,:),objk(k),conk(k,:))

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
        ymi_pos = model_u(u_pos);
        ypi_pos = plant(u_pos);
        
        % negative shift
        u_neg = uk(k-1,:) - du(i,:);
        ymi_neg = model_u(u_neg);
        ypi_neg = plant(u_neg);
        
        % get gradients
        dobjpdu(i) =   (objFun(u_pos,ypi_pos)-objFun(u_neg,ypi_neg))/(2*du(i,i));
        dconpdu(i,:) = (conFun(u_pos,ypi_pos)-conFun(u_neg,ypi_neg))/(2*du(i,i));
        
        dobjdu(i) =   (objFun(u_pos,ymi_pos)-objFun(u_neg,ymi_neg))/(2*du(i,i));
        dcondu(i,:) = (conFun(u_pos,ymi_pos)-conFun(u_neg,ymi_neg))/(2*du(i,i));
    end
    
    % modified functions
    u_last = 0;
    mod = conk(k-1,:) - conFun(uk(k-1,:),model_u(uk(k-1,:)));
    
    xObj = @(u)(optFun(u,@(u)model_u(u),objFun)+(u-uk(k-1,:))*(dobjpdu-dobjdu));
    xCon = @(u)(optFun(u,@(u)model_u(u),conFun)+(u-uk(k-1,:))*(dconpdu-dcondu)+mod);
    
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
    fprintf('%8i %10.3f %10.3f %10.3f %10.4f %10.4f %8i %10.4f\n',k,uk(k,:),objk(k),conk(k,:),flag,t1-t0)
    t0 = t1;
end
fprintf('Finished MA run [%8.2f s]\n\n',t1)

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