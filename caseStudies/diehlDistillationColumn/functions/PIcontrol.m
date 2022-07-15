% calculates the value of the controller outputs using  two individual PI
% control loops, one linking the setpoint of the temperature of tray 28 to
% the reflux return, and another linking the setpoint of the temperature of
% tray 14 to the reboiler heat duty.
function [u,de] = PIcontrol(r,K,e,y)
% ----- INPUTS --------
% r             Controller set points
% K             Controller gain
% e             Current intergral error
% y             Current value of system
% ----- OUTPUTS -------
% u             Controller output
% de            Derivative of error
% ---------------------

% Calculate control law
u = [0,0];
de = [0,0];

% Lvol controller
Kp28 = K(1);
Ki28 = K(2);
T28s = r(1);
u(1) = (4 + Kp28*(T28s-y(1)) + Ki28*e(1))/3600;

% Q controller
Kp14 = K(3);
Ki14 = K(4);
T14s = r(2);
u(2) = 2.5 + Kp14*(T14s-y(2)) + Ki14*e(2);

% Error differential
de(1) = r(1) - y(1);
de(2) = r(2) - y(2);
end