function [con] = WOconFun2(~,y)
% Calculate the constraint function for the WO case study (two-constraint 
% case used to compare PMAi and PMAj)
% --------
% y         1-by-6      outputs (yA,yB,yC,yP,yE,yG)
%
% con       1-by-1      constraint function value
% --------

con = [y(:,6) - 0.08, y(:,5)-y(:,2)+0.15];

end