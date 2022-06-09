function Q = WOquadFun(u,epsilon,c,dcdu)
% Calculates the quadratic constraint
% ---------
% u             1-nu            change in process input
% epsilon       1-nu            curavture (if 0, return nothing)
% c             1-nc            plant value
% dcdu          nu-1-nc         plant gradient
%
% Q             1-nc            quadratic value
% ---------

if epsilon == 0
    Q = [];
    return
end

Q = u*diag(epsilon)*u' + c + u*permute(dcdu,[1,3,2]);

end

