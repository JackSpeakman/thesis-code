function [obj] = WOobjFun(u,y)
% Calculate the objective function for the WO case study
% --------
% u         1-by-3      Inputs to the WO CSTR (FA,FB,TR)
% y         1-by-6      Outputs (yA,yB,yC,yP,yE,yG)
%
% obj       1-by-1      Objective function value
% --------

F = u(:,1)+u(:,2);
obj = -((1143.38.*y(:,4) + 25.92.*y(:,5)).*F - 76.23.*u(:,1) - 114.34.*u(:,2));

end

