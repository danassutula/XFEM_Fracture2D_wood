function [s,d] = LgBasis_omg(p,n)
%==========================================================================

x = p(1);
y = p(2);

switch n
    case 3
        s(1,3) = 0;
        d(2,3) = 0;
        
        s(1) = 1-x-y;
        s(2) = x;
        s(3) = y;
        
        d(1,1) =-1;
        d(1,2) = 1;
        d(1,3) = 0;
        d(2,1) =-1;
        d(2,2) = 0;
        d(2,3) = 1;
    case 6
        s(1,6) = 0;
        d(2,6) = 0;
        
        s(1) = 1-3*(x+y)+2*(x+y)^2;
        s(2) = 4*x*(1-x-y);
        s(3) = x*(2*x-1);
        s(4) = 4*x*y;
        s(5) = y*(2*y-1);
        s(6) = 4*y*(1-x-y);
        
        d(1,1) = 4*(x+y)-3;
        d(1,2) = 4*(1-2*x-y);
        d(1,3) = 4*x-1;
        d(1,4) = 4*y;
        d(1,5) = 0;
        d(1,6) =-4*y;
        d(2,1) = 4*(x+y)-3;
        d(2,2) =-4*x;
        d(2,3) = 0;
        d(2,4) = 4*x;
        d(2,5) = 4*y-1;
        d(2,6) = 4*(1-x-2*y);
end

%==========================================================================

end