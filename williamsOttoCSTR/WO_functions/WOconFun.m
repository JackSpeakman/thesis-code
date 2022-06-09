function [con] = WOconFun(~,y)
% Calculate the objective function for the WO case study
% --------
% y         1-by-6      Outputs (yA,yB,yC,yP,yE,yG)
%
% obj       1-by-1      Objective function value
% --------

con = (y(:,6) - 0.08);

end

