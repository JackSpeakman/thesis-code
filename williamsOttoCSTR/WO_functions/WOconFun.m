function [con] = WOconFun(~,y)
% Calculate the constraint function for the WO case study
% --------
% y         1-by-6      outputs (yA,yB,yC,yP,yE,yG)
%
% obj       1-by-1      constraint function value
% --------

con = y(:,6) - 0.08;

end

