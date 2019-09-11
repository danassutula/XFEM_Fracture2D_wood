function [fh_C, fh_dCdx_n] = ...
    constitutive_interpolation(C0, C1, fh_F, fh_dFdx, fh_dFdy)

if ~ismatrix(C0) || any(size(C0) ~= [3,3])
   error('`C0` must be a (3,3) matrix')
end

if ~ismatrix(C1) || any(size(C1) ~= [3,3])
   error('`C1` must be a (3,3) matrix') 
end

delta_C = C1 - C0;

function value = C(x)
    value = C0 + delta_C * fh_F(x(1), x(2));
end

function value = dCdx_n(x, n)
    value = delta_C * (fh_dFdx(x(1), x(2))*n(1) ...
                      +fh_dFdy(x(1), x(2))*n(2));
end

fh_C = @C;
fh_dCdx_n = @dCdx_n;
    
end