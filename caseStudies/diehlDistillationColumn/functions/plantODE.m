% plantODE for the distillation column case study
% Runs the distillation column with the DAE described in Dielh, M's PhD
% thesis, 2001. This is a 165 state model (166 with variable condenser
% temperature), with 82 differential equations and 83 algebraic equations.
% The liquid flow is directly calculated using the Francis wier formula.

function [dx] = plantODE(t,x,u,Fin,para)
% ----- INPUTS --------
% t             Current time (unused)
% x             Vector of states
%               [mole fractions (42), molar holdup (40), vapour flow (40), temperature (41)]
% u             Inputs [Lvol, Q];
% Fin           Feed [Fvol, x_F, T_F]
% para          Vector of parameters
% ----- OUTPUTS -------
% dx            Derivative of states
% ---------------------

%% 1. Separate inputs into variables
% parameters
nv0 = para(1);         % Reboiler holdup
nvC = para(2);         % Condenser holdup
nref = para(3);        % Francis weir parameter
eff_s = para(4);       % Stripper efficiency
eff_r = para(5);       % Rectifier efficiency
W = para(6);           % Francis weir parameter
Qloss = para(7);       % Reboiler heat loss
Ptop = para(8);        % Top pressure
dp_s = para(9);        % Stripper pressure drop
dp_r = para(10);       % Rectifier pressure drop
Ttop = para(11);       % Reflux temperature

% inputs
Lvol = u(1);
Q = u(2);

% inlet
Fvol = Fin(1);
xF = Fin(2);
Tf = Fin(3);

% column tray numbers
N = 41;
Nf = 22;

%% 2. Separate x into variables 

% Differential terms
% xi is the mole fractions of the reboiler, tray i, condenser
xi = x(1:42)';

% ni is the molar holdup of tray i [kmol]
ni = x(43:82)';

% Algebraic terms
% Vi is the vapour flow rate out of tray i [kmol/s]
Vi = [0,x(83:122)'];

% Ti is the temperatures of the reboiler, tray i, condenser [C]
Ti = [x(123:163)',Ttop];

% Fout
Fouti = x(164:165);

%% 3. Evaluate column properties
% Calculate pressures
dP = zeros(1,N+1);
dP(Nf:N) = dp_r;
dP(1:Nf-1) = dp_s;

P = ones(1,N+1)*Ptop + cumsum(dP(end:-1:1));
P = P(end:-1:1); % [Pa]

% calculate vapour concentrations [antione equation]
ant = [23.48,3626.6,-34.29;22.437,3166.4,-80.15];
Psat_f = @(T)exp(ant(:,1)-ant(:,2).*(1./(T+ant(:,3)+273.15)));

Psati = Psat_f(Ti); % [Pa]
dPsati = (ant(:,2)./(Ti+ant(:,3)+273.15).^2).*Psati; %[Pa/K]

% calculate molar volumes
coeff = [2.288/1000,0.2685,512.4,0.2453;1.235/1000,0.27136,536.4,0.24];
Vmk = 1./(coeff(:,1)).*coeff(:,2).^(1+(1-(273.15+Ti)./coeff(:,3)).^coeff(:,4)); %[L/kmol]
Vm = dot([xi;(1-xi)],Vmk); % [L/kmol]

nvi = ni.*Vm(2:N); % volume holdup [L]
n0 = nv0/Vm(1);    % molar hold up [kmol]
nC = nvC/Vm(N+1);  % molar hold up [kmol]

% calculate liquid flowrate
Li = [0, W*(nvi-nref).^(3/2)];  % [L/s]
Li = real(Li)./Vm(1:N);         % [kmol/s]

%Vi(Vi<0) = 0;
%Li(Li<0) = 0;

Vmkf =  1./(coeff(:,1)).*coeff(:,2).^(1+(1-(273.15+Tf)./coeff(:,3)).^coeff(:,4)); %Feed
Vmf = dot([xF;(1-xF)],Vmkf); %[L/kmol]

a = (1-(Ti+273.15)./coeff(:,3));
dVmkdT = -(coeff(:,4)./(coeff(:,1).*coeff(:,3))).*log(coeff(:,2)).*a.^(coeff(:,4)-1).*coeff(:,2).^(1+a.^(coeff(:,4)));
dVmdT = dot([xi;(1-xi)],dVmkdT);    % [kmol/L/K]
dVmdX = Vmk(1,:) - Vmk(2,:);        % [kmol/L]

% Calculate vapour compositions
eff = zeros(1,N);
eff(Nf+1:N) = eff_r;
eff(2:Nf) = eff_s;
eff(1) = 1;
yi = zeros(1,N);
yi(1) = Psati(1,1)/P(1,1)*x(1,1);
for i = 2:N
    yi(1,i) = eff(1,i).*(Psati(1,i)./P(1,i)).*xi(1,i) + (1-eff(1,i))*yi(1,i-1);
end

% calculate enthalpies
hki = [18.31,1.713e-2,6.399e-5;
    31.92,4.49e-2,9.663e-5];
Tc = [512.6;536.7];
Pc = [8.096e6;5.166e6];
Omega = [0.557;0.612];
a = 6.09648;
b = 1.28862;
c = 1.016;
d = 15.6875;
e = 13.4721;
f = 2.615;

Tr = (1./Tc)*(Ti+273.15);
Pr = (1./Pc)*P;

hlk = 4.186*(hki(:,1)*(Ti)+hki(:,2)*(Ti.^2)+hki(:,3)*(Ti.^3));
hvk = hlk + 8.314*Tc.*sqrt(1-Pr.*Tr.^-3).*(a-b*Tr+c*Tr.^7+Omega.*(d-e*Tr+f*Tr.^7));
hfk = 4.186*(hki(:,1)*(Tf)+hki(:,2)*(Tf.^2)+hki(:,3)*(Tf.^3));

hl = dot([xi;(1-xi)],hlk);          % [kJ/kmol]
hv = dot([yi;(1-yi)],hvk(:,1:N));   % [kJ/kmol]
hf = dot([xF;(1-xF)],hfk);          % [kJ/kmol]

dhldX = hlk(1,:) - hlk(2,:);
dhlkdT = 4.186*(hki(:,1)+hki(:,2)*Ti*2+hki(:,3)*(Ti.^2)*3);
dhldT = dot([xi;(1-xi)],dhlkdT);


%% ###################### SOLVE BALANCES ######################
x_dot = zeros(size(xi));

%% 4. Solve reboiler balance to calculate V0
c1 = -(Psati(1,:)-Psati(2,:))./(dPsati(1,:).*xi+dPsati(2,:).*(1-xi));
c2 = (dVmdX(1)+c1(1)*dVmdT(1));
c3 = -n0*c2/Vm(1);
c4 = [n0,ni,nC].*(dhldX+dhldT.*c1);

% linear equation [dx, dn, V, B];
A = [-c3,1,0,0;                     % Constant volume
    0,1,1,1;                        % Overall mass conservation
    n0,xi(1),yi(1),xi(1);           % Component mass conservation
    c4(1),hl(1),hv(1),hl(1)];       % Energy balance

B = [0;
    Li(2);
    Li(2)*x(2);
    (Q-Qloss)+Li(2)*hl(2)];

X = A\B;

x_dot(1) = X(1);
n_dot0 = X(2);
V0 = X(3);
Bott = X(4);

Vi(1) = V0;

%% 5. Solve condenser balance
c2 = (dVmdX(N+1)+c1(N+1)*dVmdT(N+1));
c3 = -nC*c2/Vm(N+1);
Lcond = Lvol/Vm(N+1);           % [kmol/s]

% linear equation [dx, dn, D]
A = [-c3, 1, 0;
    0, 1, 1;
    nC, xi(N+1), xi(N+1)];

B = [0;
    Vi(N) - Lcond;
    Vi(N)*yi(N) - Lcond*xi(N+1)];

X = A\B;

x_dot(N+1) = X(1);
n_dotC = X(2);
Dist = X(3);

%% 6. Solve molar hold up balance
Li(N+1) = Lcond;
Fi = zeros(1,N);
Fi(Nf) = Fvol/Vmf;

n_dot = Vi(1:N-1)+Li(3:N+1)-Vi(2:N)-Li(2:N)+Fi(2:N);        % Trays

%% 7. Solve component balance

x_dot(2:N) = (Vi(1:N-1).*yi(1:N-1) - Vi(2:N).*yi(2:N) + Li(3:N+1).*xi(3:N+1) - (Li(2:N)+n_dot).*xi(2:N) + Fi(2:N).*xF)./ni;

%% 8. Solve Algebraic equations
%L_dot = W*(nvi-nref).^(3/2);

V_dot = (1./ni).*(Vi(1:N-1).*hv(1:N-1) - Vi(2:N).*hv(2:N) + Li(3:N+1).*hl(3:N+1) - (Li(2:N)+n_dot).*hl(2:N) + Fi(2:N).*hf) - ...
    c4(2:N).*x_dot(2:N);
T_dot = (P(1:N) - Psati(1,1:N).*xi(1:N) - Psati(2,1:N).*(1-xi(1:N)))./P(1:N);
Fout_dot = [Dist-Fouti(1), Bott-Fouti(2)];

dx = [x_dot,n_dot,V_dot,T_dot,Fout_dot]';

if ~isreal(dx)
    warning('not real')
end

end

