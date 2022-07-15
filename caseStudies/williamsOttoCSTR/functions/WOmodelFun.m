function [ySol] = WOmodelFun(u,yGuess,dtheta)
% WOmodelFun takes the inputs to the WO model and calculate the outputs
% --------
% u         1-by-3      Inputs to the WO CSTR (FA,FB,TR)
% yGuess    1-by-6      Guess of outputs
% dtheta    1-by-2      Uncertainty in parameters
%
% ySol      1-by-6      Outputs (yA,yB,yC,yP,yE,yG)
% --------

if nargin < 2 || isempty(yGuess) % no guess
    yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
end

% inputs
F_Ain = u(1);
F_Bin = u(2);
T     = u(3);

% nominal model parameters
k_0 = [1.06e+07    3.1e+11]; %[2.189*1e8, 4.310*1e13]; % 1/s
E = [7.034e+03    1.068e+04]; %[8077.6,12438.5]; % K
M = 2105; % kg

if nargin == 3 % add uncertainty
    dtheta = [dtheta(:);zeros(5,1)];
    E   = E + dtheta(1:2)';
    k_0 = k_0 + dtheta(3:4)';
    M = M + dtheta(5);
end

% other variables
F  = F_Ain + F_Bin; %kg/s
k1 = k_0(1)*exp(-E(1)/(T+273.15)); %1/s
k2 = k_0(2)*exp(-E(2)/(T+273.15)); %1/s

%% mass balances
optionX = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'StepTolerance',1e-12,'FunctionTolerance',1e-12);
[ySol,~,flag] = fsolve(@(y)WOmodelODE(y,k1,k2,M,F_Ain,F_Bin,F), yGuess, optionX);

if any(ySol<0 | ySol>1) % bad solve - try new yGuess
    ySol = fsolve(@(y)WOmodelODE(y,k1,k2,M,F_Ain,F_Bin,F), [0.5,0.5,0,0,0,0], optionX);
end
    
if flag ~= 1
    % pause if issue
    keyboard
end

end

