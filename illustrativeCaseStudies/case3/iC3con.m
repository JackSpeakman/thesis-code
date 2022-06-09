function con = iC3con(u1,u2,~,th2,i)
%iC3con is the constraint function for the third illustrative case.
% ----------
% u1    n-by-1      input 1
% u2    n-by-1      input 2 
% th1   n-by-1      parameter 1 [unused]
% th2   n-by-1      parameter 2
% i     [0,1,2]     model request
%
% con   n-by-1      constraint function value
% ----------

% calculate zeta
if i == 0 % plant
    zeta = exp(u1)+exp(-u1)-1;
elseif i == 1 % model 1
    zeta = u1.^2+1;
elseif i == 2 % model 2
    zeta = sqrt(3*u1.^2+1);
end

% calculate constraint
con = th2*(zeta-0.9*u2).^2+u2-1.8;

end