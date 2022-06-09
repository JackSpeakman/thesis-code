function con = iC2con(u1,u2,th1,th2,i)
%iC2con is the constraint function for the second illustrative case.
% ----------
% u1    n-by-1      input 1
% u2    n-by-1      input 2 
% th1   n-by-1      parameter 1
% th2   n-by-1      parameter 2
% i     [~,1,2]     indicator for the constraint
%
% con   n-by-2      constraint
% ----------

% calculate the constraint

if ~exist('i','var')
    con1 = zeros([size(u1)]);
    con2 = zeros([size(u1)]);
    con1 = th1.*(u1).^2+th1.*(u2).^2-5*u1+u2+0.5;
    con2 = th2.*(u1).^2+th2.*(u2).^2+1.2*u1+1.2*u2-3;
    con = cat(sum(size(u1)>1)+1,con1,con2);
elseif i == 1
    con = th1.*(u1).^2+th1.*(u2).^2-5*u1+u2+0.5;
else
    con = th2.*(u1).^2+th2.*(u2).^2+1.2*u1+1.2*u2-3;
end    

end