function [con] = iC1con(u,th)
%illCase1con is the constraint function for the first illustrative case.
% This case study is a one-dimensional problem with a single input variable
% of (u), and a single constraint (no objective).
% ----------
% u     1-by-n      inputs to constraint
% th    m-by-2      parameters to constraint
%
% con   n-by-m      constraint
% ----------

% correct input arguments
u = u(:)';
if size(th,2)~=2
    if size(th,1)~=2
        error('Size of theta incorrect')
    end
    th = th';
end

% calculate the constraint
con = (u-2).^2-3-th(:,1).*sin(th(:,2)*u);

end

