function [model,plant,objFun,conFun] = setupDCshortFunctions
% Set-up variables to be given to run__MA. The decision variables are the 
% controller setpoints (r). four functions are required (model(r), 
% plant(r), objFun(r,y), conFun(r,y)). The cost/cons functions are defined 
% in terms of the input variables (u). PIcontoller relates r/y to u. The
% following code sets up the distillation column for all models (variable
% size of rectifier, but only uses the nominal (Nr = 9) for the standard
% MA run.
% ------ OUTPUT VARIABLES ------
% model     @(r,i)          Model function.
%                               r       decision variable
%                               i       DC variable parameter {1,2,3,4,5}
% plant     @(r)            Plant function.
%                               r       decision variable
% objFun    @(r,y)          Objective function.
%                               r       decision variable
%                               y       DC ODE output
% conFun    @(r,y)          Constraint function.
%                               r       decision variable
%                               y       DC ODE output
%
% ------ EXAMPLES ------
% >> [model,plant,objFun,conFun] = setupDCshortFunctions;
%       This will create the functions required for running the MA schemes
%
% See '../../../methods/modifierAdaptation/distillationColumn/makePlantsMA_DCshort.m'
%   for a script using this function

%% 0. General set-up
% define feed
Fin = [14/3600,0.32,71];
Fdiff = [32/30,30/32,1];

% time function used by ODE solver
tFun = @(t)(logspace(2,log10(t+100),10001)-100);
tspan = 60000000;

% controller kTT is the keyTrayTemp function (defined later)
K = [-0.1,-0.0002,0.1,0.0002];                         % controller gain
cont = @(r,y,kTT)(PIcontrol(r,K,y,kTT(y)).*[3600,1]);  % returns u

% ODE output rearranger
% plant/model outputs end in [F_bott, F_dist, x_bott, x_dist] which allows
% for the length of x to vary without affecting the solution of the obj/con
% function evaluation.
ODEout = @(r,sol,Nout,kTT)([sol.y(:,end);...    % standard outputs
    cont(r,sol.y(:,end),kTT)';...               % inputs (Lvol,Q)
    sol.y(end-1,end);sol.y(end,end);...         % dist/bott flowrate
    sol.y(Nout(1),end);sol.y(Nout(2),end)]);    % dist/bott conc.


%% 1. Plant function
% parameters to be given to 'plantPI' and 'plantODE'
theta_p = [8.5,...      % Reboiler holdup
    0.17,...            % Condenser holdup
    0.155,...           % Francis weir parameter
    0.62,...            % Stripper efficiency
    0.35,...            % Rectifier efficiency
    0.166,...           % Francis weir parameter
    0.5,...             % Reboiler heat loss
    93900,...           % Top pressure
    250,...             % Stripper pressure drop
    190,...             % Rectifier pressure drop
    47.2];              % Reflux temperature

% key tray temperature function (extracts from plantPI output)
keyTrayTempp = @(y)(y(125+[28,14]));

% ODE solver set-up
M = diag([ones(1,84),zeros(1,83)]); % ODE mass matrix
opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-8);
plantsol_y0 = @(r,y0)ode15s(@(t,y)plantPI(t,y,r,Fin.*Fdiff,K,theta_p),tFun(tspan),y0,opts);

% get outputs from ODE solution
Noutp = [3,44];         % x index of bottoms/distillate conc.
plant_y0 = @(r,y0)(ODEout(r,plantsol_y0(r,y0),Noutp,keyTrayTempp));

% remove x0 from function by running plant at some point
x0 = [linspace(0.001,0.01,5),linspace(0.012,0.988,32),linspace(0.99,0.999,5)];
n0 = linspace(0.0034,0.0044,40);
V0 = linspace(0.0028,0.0032,40);
T0 = linspace(95,65,41);
Fout0 = [4/3600,10/3600]/69;

r0 = [71.7,84];
yp0 = plant_y0(r0,[0,0,x0,n0,V0,T0,Fout0]');

plant = @(r)plant_y0(r,yp0(1:end-6));

DCSolver(@(t,y)plantPI(t,y,r0,Fin.*Fdiff,K,theta_p),yp0(1:end-6),tFun(tspan),opts,@(y)(y(end-1:end)),@(y)(y(Noutp)),keyTrayTempp,@(y,T)(PIcontrol(r0,K,y,T).*[3600,1]));

%% 2. Model function
% Define parameters (i \in {1,2,3,4,5} is the surrogate theta for the DC)
Ns = 14;                % number of stripper trays
Nr = 7:11;              % number of rectifier trays

N = @(i)(Ns+Nr(i));     % total number of trays
Nf = Ns;                % feed tray
mu = 3;                 % nominal model number

% parameters to be given to 'modelPI_short' and 'modelODE_short'
theta_m = @(i)[8.5,...  % Reboiler holdup
    0.17,...            % Condenser holdup
    0.155,...           % Francis weir parameter
    1,...               % Stripper efficiency [plant - 0.62]
    1,...               % Rectifier efficiency [plant - 0.35]
    0.166,...           % Francis weir parameter
    0.51,...            % Reboiler heat loss [plant - 0.5]
    90900,...           % Top pressure [plant - 93900]
    250,...             % Stripper pressure drop
    190,...             % Rectifier pressure drop
    44.2,...            % Reflux temperature [plant - 47.2]
    N(i),...            % Number of trays
    Nf];                % Feed tray

% key tray temperature function (extracts from plantPI output)
keyTrayTempm = @(y,i)(y(2*N(i)+2+[N(i)-4,7]));

% ODE solver set-up
M = @(i)(diag([ones(1,N(i)+3),zeros(1,2*N(i)+3)])); % mass matrix
opts = @(i)(odeset('Mass',M(i),'RelTol',1e-10,'AbsTol',1e-10));
modelsol = @(r,y0,i)ode15s(@(t,y)modelPI_short(t,y,r,Fin,K,theta_m(i)),tFun(tspan),y0,opts(i));

% get outputs from ODE solution (with y0)
Noutm = @(i)([3,N(i)+3]);
model_y0 = @(r,y0,i)(ODEout(r,modelsol(r,y0,i),Noutm(i),@(y)(keyTrayTempm(y,i))));

% get the initial conditions for each model (as size of y changes with Nr)
y0 = cell(numel(Nr),1);

for ii = 1:numel(Nr)
    % run each model at random initial point
    x0 = linspace(0.01,0.99,N(ii)+1);
    V0 = linspace(0.0028,0.0032,N(ii));
    T0 = linspace(95,65,N(ii)+1);
    Fout0 = [4/3600,10/3600]/69;
    
    r0 = [71.7,84];
    y0{ii} = model_y0(r0,[0,0,x0,V0,T0,Fout0]',ii);
end

% finally get the model from the basic inputs (without y0)
model = @(r,i)(model_y0(r,y0{i}(1:end-6),i));

%% 3. Objective function
% set-up function for cost (w.r.t. u and y)
Profit = [200*32,400*60]*3.6;

% Q is heat flux, F_Bott is bottoms flowrate, F_Dist is distilate flowrate
cost = @(Q,F_Dist,F_Bott)(Q*0.12-F_Dist*Profit(1)-F_Bott*Profit(2)+3.7);

% obj in terms of r and y
objFun = @(r,y)(cost(y(end-4),y(end-3),y(end-2)));

%% 4. Constraint function
% Set-up function for constraint (w.r.t. u and y)
cons = @(x_Bott,x_Dist)([x_Bott-0.005,0.986-x_Dist]);

% con in terms of r and y
conFun = @(r,y)(cons(y(end-1),y(end)));

end








