function [dF] = WOmodelODE(y,k1,k2,M,F_Ain,F_Bin,F)

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