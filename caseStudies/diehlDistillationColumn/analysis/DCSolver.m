function [y,t,yt] = DCSolver(ODE,y0,tspan,opts,Fout,Nout,kTT,cont)
% DCSolver takes the ODE for a DC (either plantPI or modelPI_short) and
% solves it using ode15s, returning the rearranged output which can be 
% given to the objFun and conFun no matter the ODE given.
% 
% ------ INPUT VARIABLES ------
% r         1-by-2          Controller setpoints
% ODE       @(t,y)          ODE function to be solved
% y0        1-by-n_y        Initial conditions
% tspan     1-by-n_t        Time span for ODE solver
% opts      struct          odeset options
% Fout      @(y)            Gets the dist/bott flowrates from y
% Nout      @(y)            Gets the dist/bott concentrations from y
% kTT       @(y)            Gets the key tray temperatures from y
% cont      @(y,T)          Solves the controller
%
% ------ OUTPUT VARIABLES ------
% y     1-by-(n_y+6)    Outputs of ODE
%                           (plus controller solution (i.e. u)
%                           dist/bott flowrates,
%                           dist/bott conc.)
% t     1-by-n_t2       Time steps
% yt    n_y-by-n_t2     Outputs of ODE at each time step (for plotting)
% 
% ------ EXAMPLES ------
% See 'setupDCshortFunctions.m' examples

%% Run ODE
sol = ode15s(ODE,tspan,y0,opts);

% solve outputs
t = sol.x;
yt = sol.y;

%% Get extended outputs
% get dist/bott flowrates from y
F = Fout(yt(:,end));

% get dist/bott concentrations from y
N = Nout(yt(:,end));

% get key tray temperatures from y
T = kTT(yt(:,end));

% get inputs
u = cont(yt(:,end),T)';

% combine
y = [yt(:,end); u; F; N];

end