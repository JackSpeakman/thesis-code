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
[ySol,~,flag] = fsolve(@(y)massBalance(y), yGuess, optionX);

if any(ySol<0 | ySol>1) % bad solve - try new yGuess
    ySol = fsolve(@(y)massBalance(y), [0.5,0.5,0,0,0,0], optionX);
end
    
if flag ~= 1
    a=1;
end


    function [dF] = massBalance(y)
        
        R1 = k1*y(1)*y(2);
        R2 = k2*y(2)*y(3);
        R3 = k3*y(3)*y(4);
        %       IN    - OUT    +GEN/-CON R1   +GEN/-CON R2
        dF = [  F_Ain - F*y(1) -   M*R1                    ; %A
                F_Bin - F*y(2) -   M*R1 -   M*R2           ; %B
                0     - F*y(3) + 2*M*R1 - 2*M*R2 -     M*R3; %C
                0     - F*y(4)          +   M*R2 - 0.5*M*R3; %P
                0     - F*y(5)          + 2*M*R2           ; %E
                0     - F*y(6)                   + 1.5*M*R3]; %G
    end

end