function [ySol] = WOplantFun(u,yGuess,dtheta)
% WOmodelFun takes the inputs to the WO model and calculate the outputs
% --------
% u         1-by-3      Inputs to the WO CSTR (FA,FB,TR)
% yGuess    1-by-6      Guess of outputs
% dtheta    1-by-3      Uncertainty in parameters
%
% ySol      1-by-6      Outputs (yA,yB,yC,yP,yE,yG)
% --------
if nargin < 2 || isempty(yGuess) % no guess
    yGuess = [0.141, 0.332, 0.023, 0.103, 0.284, 0.117];
end

% inputs
F_Ain = u(1);
F_Bin = u(2);
T     = u(3);

% nominal model parameters
k_0 = [1.6599e6, 7.2117e8, 2.6745e12]; %1/s
E = [6666.7, 8333.3, 11111]; %1/K
M = 2105; % kg
if nargin > 2
    E = E + dtheta;
end

% other variables
F  = F_Ain + F_Bin; %kg/s
k1 = k_0(1)*exp(-E(1)/(T+273.15)); %1/s
k2 = k_0(2)*exp(-E(2)/(T+273.15)); %1/s
k3 = k_0(3)*exp(-E(3)/(T+273.15)); %1/s

%% mass balances
optionX = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'StepTolerance',1e-12,'FunctionTolerance',1e-12);
[ySol,~,flag] = fsolve(@(y)WOplantODE(y,k1,k2,k3,M,F_Ain,F_Bin,F), yGuess, optionX);

if any(ySol<0 | ySol>1) % bad solve - try new yGuess
    ySol = fsolve(@(y)WOplantODE(y,k1,k2,k3,M,F_Ain,F_Bin,F), [0.5,0.5,0,0,0,0], optionX);
end
    
if flag ~= 1
    keyboard;
end




end