function obj = iC3obj(u1,u2,th1,~,i)
%iC3obj is the objective function for the third illustrative case.
% ----------
% u1    n-by-1      input 1
% u2    n-by-1      input 2 
% th1   n-by-1      parameter 1
% th2   n-by-1      parameter 2 [unused]
% i     [0,1,2]     model request
%
% obj   n-by-1      objective function value
% ----------

% calculate zeta
if i == 0
    zeta = exp(u1)+exp(-u1)-1;
elseif i == 1
    zeta = u1.^2+1;
elseif i == 2
    zeta = sqrt(3*u1.^2+1);
end

% calculate objective
obj = 4-th1*zeta - u2;

end