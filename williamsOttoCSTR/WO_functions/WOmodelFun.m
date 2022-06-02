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
[ySol,~,flag] = fsolve(@(y)massBalance(y), yGuess, optionX);

if any(ySol<0 | ySol>1) % bad solve - try new yGuess
    ySol = fsolve(@(y)massBalance(y), [0.5,0.5,0,0,0,0], optionX);
end
    
if flag ~= 1
    a=1;
end

    function [dF] = massBalance(y)
        
        R1 = k1*y(1)*y(2)^2;
        R2 = k2*y(1)*y(2)*y(4);
        %       IN    - OUT    +GEN/-CON R1   +GEN/-CON R2
        dF = [  F_Ain - F*y(1) -   M*R1       -   M*R2; %A
                F_Bin - F*y(2) - 2*M*R1       -   M*R2; %B
                0     - F*y(3)                        ; %C
                - F*y(4) +   M*R1       -   M*R2; %P
                - F*y(5) + 2*M*R1               ; %E
                - F*y(6)                + 3*M*R2]; %G
        
    end

end

