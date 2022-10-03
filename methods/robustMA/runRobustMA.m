function [uk,yk,conk,objk] = runRobustMA(varargin)
% The master code for running all the robust modifier adaptation methods.
% The robust MA approaches are described in the I&EC article titled "Robust 
% Modifier Adaptation via Worst-Case and Probabilistic Approaches" by Jack 
% Speakman and Gregory Francois. This framework describes three different
% methods of improving the robustness of the MA approach.
% 
% ------ INPUT VARIABLES ------
% varargin                                  cell of inputs
%   'filter'            1-by-1              Input filter
%   'kmax'              1-by-1              Number of RTO iterations
%   'stratingPoint'     1-by-n_u            RTO starting point
%   'uGuess'            1-by-n_u            Input guess (for initial runs)
%   'th'                i-by-n_th           Set of model parameters
%   'th_nom'            i-by-n_th           Nominal model parameters
%   'conFun'            @(u,y)              Constraint function
%   'objFun'            @(u,y)              Objective function
%   'modelFun'          @(u,th)             Model function
%   'method'            string              ['WCMA','PMAi','PMAj']
%   'umin'              1-by-n_u            Optimization minimum limit
%   'umax'              1-by-n_u            Optimization maximum limit
%   'probChance'        1-by-1              Probability constraint chance
% 
% ------ OUTPUT VARIABLES ------
% uk        kmax-by-n_u         Inputs for iterations 1 to kmax
% yk        kmax-by-n_y         Outputs for iterations 1 to kmax
% conk      kmax-by-n_c         Constraints for iterations 1 to kmax
% objk      kmax-by-1           Objective for iterations 1 to kmax
% 
% ------ EXAMPLES ------
% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runRobustMA
%       This runs the default argument which is WCMA applied to WO
%
% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runRobustMA('method','PMAi')
%       This runs the same set-up as the default but uses the PMAi method
% 
% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runRobustMA('method','PMAj','filter',0.7,'startingPoint',[3.2,6.5,90])
%       This runs the PMAj method to WO with an input filter gain of 0.7, 
%       and the initial point of [3.2,6.5,90]
% 
% See 'distillationColumn/makePlotsWCMA_DC.m' for distillation column example
% See 'semiBatchReactor/makePlotsWCMA_DC.m' for semi-batch reactor example


%% 0. Deal with varargin
% set-up default values [for WCMA on WO]
% general RTO parameters
filter = 0.5;                       % input filter
kmax = 11;                          % number of iterations
u0 = [];                            % starting point

% model/plant parameters
th_nom = [0,0];                     % default nominal parameters
th = [0,0; 70,160; -70,160; 70,-160; -70,-160]*3; % set of model parameters
conFun = @(u,y)WOconFun(u,y);       % constraint function
objFun = @(u,y)WOobjFun(u,y);       % objective function
model_th = @(u,th)WOmodelFun(u,[],th); % model function [with paramter input]
plant = @(u)WOplantFun(u);          % plant function
umin = [3,6,80];                    % optimization limit [min]
umax = [4.5,11,105];                % optimization limit [max]

% method parameters
method = 'WCMA';                    % robust MA method ['WCMA','PMAi','PMAj']
rho = 0.95;                         % probability constraint chance

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
            objFun = varargin{i+1};
        case 'modelFun'
            model_th = varargin{i+1};
        case 'plantFun'
            plant = varargin{i+1};
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
fprintf('\n#### %s - %4.2f ####\n\n',method,filter);

%% 1. Set-up parameters
tic

% uGuess
if ~exist('uGuess','var')
    uGuess = (umin+umax)/2;
end
n_u = numel(uGuess);

% combined model function
yGuess = model_th(uGuess,th_nom);

n_th = size(th,1);
th = [th_nom;th];

model = @(u,i)model_th(u,th(i+1,:));

% fix conFun size
size_c = size(conFun(uGuess,yGuess));
n_c = prod(size_c);
if size_c(1) ~= 1
    conFun = @(u,y)(conFun(u,y)');
end

% get model outputs size
n_ym = numel(yGuess);

yGuessp = plant(uGuess);
n_yp = numel(yGuessp);

% gradient shift
du = eye(n_u)*0.0001;

% proallocate outputs
uk = zeros(kmax,n_u);
yk = zeros(kmax,n_yp);
objk = zeros(kmax,1);
conk = zeros(kmax,n_c);

% global variables
u_last = zeros(1,n_u);
y_last = cell(1,n_th);

%% 2. Find starting point
k = 1;
fminconopts = optimoptions('fmincon','Display','off');
fminimaxopts = optimoptions('fminimax','Display','off','MaxFunctionEvaluations',10000,'MaxIterations',1000);

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
fprintf('Beginning %s run\n',method)
fprintf(['%8s ' repmat('%10s ',1,n_u) '%10s ' repmat('%10s ',1,n_c) '%8s %10s [s]\n'],...
    'k',string([repmat('Input ',n_u,1),num2str((1:n_u)')]),'Cost',string([repmat('Con ',n_c,1),num2str((1:n_c)')]),'Flag','Time');

fprintf(['%8s ' repmat('%10.3f ',1,n_u) '%10.4f ' repmat('%10.3f ',1,n_c) '\n'],'init',uk(k,:),objk(k),conk(k,:))

% run robust MA
for k = 2:kmax
    
    % gradients
    dobjpdu = zeros(n_u,1);
    dconpdu = zeros(n_u,1,n_c);
    dobjdu = zeros(n_u,n_th);
    dcondu = zeros(n_u,n_th,n_c);
    
    for i = 1:n_u
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
    
    % modifiers
    modO = zeros(n_th,1);
    modC = zeros(n_th,n_c);
    for j = 1:n_th
        modO(j,:) = objk(k-1,:) - objFun(uk(k-1,:),model(uk(k-1,:),j));
        modC(j,:) = conk(k-1,:) - conFun(uk(k-1,:),model(uk(k-1,:),j));
    end
    
    % modified functions [WCMA]
    modObjFun = @(u)(WCoptFun(u,model,objFun,1)+modO+(dobjpdu-dobjdu)'*(u-uk(k-1,:))');
    modConFun = @(u)(WCoptFun(u,model,conFun,n_c)+modC+...
        permute(sum(repmat((u-uk(k-1,:))',1,n_th,n_c).*(dconpdu-dcondu),1),[2,3,1]));
    
    switch method
        case 'WCMA'
            % objective for WCMA is expected value (mean) but can also be
            % the max. This code can be modified to deal with the max by
            % solving the 'new optimum' within the individual cases and
            % using fminimax for WCMA (this avoid issues with continuity)
            xObj = @(u)mean(modObjFun(u));
            
            % constraint for WCMA is all the constraints
            xCon = @(u)(modConFun(u));
            
        case 'PMAi'
            % objective for PMAi is expected value (mean) but can also be
            % the probability constrained objective using rho
            xObj = @(u)individualProbCalc(modObjFun(u),0.5);
            
            % constraint for PMAi is the probCalc using rho
            xCon = @(u)individualProbCalc(modConFun(u),rho);
            
        case 'PMAj'
            % objective for PMAj is expected value (mean) but can also be
            % the probability constrained objective using rho
            xObj = @(u)individualProbCalc(modObjFun(u),0.5);
            
            % constraint for PMAj probCalc using the joint approach
            xCon = @(u)jointProbCalc(modConFun(u),rho);
            
    end
    
    % new optimum
    [uOpt,~,~,flag] = fminimax(@(u)xObj(u),uGuess,[],[],[],[],umin,umax,...
        @(u)deal(xCon(u),[]),fminimaxopts);
            
    % apply filter
    uk(k,:) = uk(k-1,:)*(1-filter) + uOpt*filter;
    
    % run plant
    yk(k,:) = plant(uk(k,:));
    objk(k,:) = objFun(uk(k,:),yk(k,:));
    conk(k,:) = conFun(uk(k,:),yk(k,:));
    
    % print result
    t1 = toc;
    fprintf(['%8i ' repmat('%10.3f ',1,n_u) '%10.4f ' repmat('%10.3f ',1,n_c) '%8i %10.4f\n'],k,uk(k,:),objk(k),conk(k,:),flag,t1-t0)
    t0 = t1;
end
fprintf('Finished %s run [%8.2f s]\n\n',method,t1)

%% 4. Embedded functions
    function out = WCoptFun(u,yFun,outFun,n_out)
        % calculates the optimization function for all th (efficient)
        if u == u_last
            y = y_last;
        else
            y = cell(1,n_th);
            for jj = 1:n_th
                y{jj} = yFun(u,jj);
            end
            u_last = u;
            y_last = y;
        end
        
        out = zeros(n_th,n_out);
        for jj = 1:n_th
            out(jj,:) = outFun(u,y{jj});
        end
    end

    function pOut = individualProbCalc(out,rho)
        % calculates the probabilistic value (individual)
        mu = mean(out);
        Sigma2 = var(out);
        zp = sqrt(2)*(erfinv(2*rho-1));
        pOut = mu+sqrt(Sigma2)*zp;
    end

    function pOut = jointProbCalc(out,rho)
        % calculates the probabilistic value (joint)
        mu = mean(out);
        Sigma2 = cov(out);
        zp = sqrt(2)*(erfinv(2*rho-1));
        pOut = mu+sqrt(diag(Sigma2)')*zp;
        
        if numel(mu) > 1
            opts = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
            pOut = fsolve(@(x)(mvncdf([-inf,-inf],[0,0],mu-x,Sigma2)-rho),mu,opts);
        end
    end
end