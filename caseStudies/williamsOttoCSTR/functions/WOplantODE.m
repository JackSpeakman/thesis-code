function [dF] = WOplantODE(y,k1,k2,k3,M,F_Ain,F_Bin,F)

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