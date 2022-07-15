% Solves the distillation column for the plant using a PI controller
function [dxPI] = plantPI(t,xPI,r,Fin,K,para)
% ----- INPUTS --------
% t             Current time (unused)
% xPI           Vector of states
%               [controller intergral term (2), mole fractions (42), molar holdup (40), vapour flow (40), temperature (41)]
% r             Set points [T28, T14]
% Fin           Feed [Fvol, x_F, T_F]
% K             Controller gain
% para          Vector of parameters
% ----- OUTPUTS -------
% dxPI          Derivative of states
% ---------------------

%% 1. Initialise parameters/inputs/etc.
if ~exist('para','var') % set to nominal values (plant)
    % theta is the parameters of the model
    % [nv0,nvC,nref,eff_s,eff_r,Wtray,Qloss,Ptop,dp_s,dp_r,Tc]
    para = [8.5, 0.17, 0.155, 0.62, 0.35, 0.166, 0.51, 93900, 250, 190, 47.2];
end

% Initial states xPI = [e(2),xi(42),ni(40),Vi(40),Ti(41)]
e = xPI([1,2]);      % PI controller integral error on tray 28 and 14 (key trays)

N_key = [28,14];            % Key trays
T_key = xPI(N_key+125);      % Temperature onkey trays

% Column
N = 41;                 % Trays including reboiler
Nf = 22;                % Feed tray (+ reboiler)

%% 2. Calculate controller
[u,de] = PIcontrol(r,K,e,T_key);
dxPI = zeros(size(xPI));
dxPI(1:2) = de;

% Error differential
if u(1)<1/3600
    u(1) = 1/3600;
    dx(1) = 0;
elseif u(1)>5/3600
    u(1) = 5/3600;
    dx(1) = 0;
else
    dx(1) =r(1) - T_key(1);
end

if u(2)<1
    u(2) = 1;
    dx(2) = 0;
elseif u(2)>3
    u(2) = 3;
    dx(2) = 0;
else
    dx(2) = r(1) - T_key(2);
end

%% Run open-loop plant
dxPI(3:167) = plantODE(t,xPI(3:167),u,Fin,para);

end