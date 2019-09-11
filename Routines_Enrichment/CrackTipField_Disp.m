function [u,dudx,dudy] = CrackTipField_Disp(K1,K2,G,k,r,theta)

n  = length(r);
u  = zeros(n,2);

dudr = zeros(n,2);
dudt = zeros(n,2);
dudx = zeros(n,2);
dudy = zeros(n,2);

drdx = cos(theta);
drdy = sin(theta);

dtdx = -drdy./r;
dtdy =  drdx./r;

theta = theta/2;
c = cos(theta);
s = sin(theta);
sqrt_r = sqrt(r);

if K1 ~= 0
    
    const = (K1/2/G/sqrt(2*pi));
    
    u(:,1) = c.*(k-1+s.^2*2).*sqrt_r*const;
    u(:,2) = s.*(k+1-c.^2*2).*sqrt_r*const;
    
    const = const*0.5;
    
    dudr(:,1) = c.*(k-1+s.^2*2)./sqrt_r*const;
    dudr(:,2) = s.*(k+1-c.^2*2)./sqrt_r*const;
    
    dudt(:,1) = s.*(c.^2*6-k-1).*sqrt_r*const;
    dudt(:,2) = c.*(k+5-c.^2*6).*sqrt_r*const;
    
end

if K2 ~= 0
    
    const = (K2/2/G/sqrt(2*pi));
    
    u(:,1) = u(:,1) + s.*(k+1+c.^2*2).*sqrt_r*const;
    u(:,2) = u(:,2) - c.*(k-1-s.^2*2).*sqrt_r*const;
    
    const = const*0.5;
    
    dudr(:,1) = dudr(:,1) + s.*(k+1+c.^2*2)./sqrt_r*const;
    dudr(:,2) = dudr(:,2) + c.*(1-k+s.^2*2)./sqrt_r*const;
    
    dudt(:,1) = dudt(:,1) + c.*(k-3+c.^2*6).*sqrt_r*const;
    dudt(:,2) = dudt(:,2) + s.*(k-3+c.^2*6).*sqrt_r*const;
    
end

for i = 1:n
   
    dudx(i,:) = dudr(i,:)*drdx(i) + dudt(i,:)*dtdx(i);
    dudy(i,:) = dudr(i,:)*drdy(i) + dudt(i,:)*dtdy(i);
    
end
