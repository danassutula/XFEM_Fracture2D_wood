function [D, alpha, kappa] = MtxConstit(n, E, v, probtype)

D = cell(n,1);

alpha = zeros(n,1); % Fracture mechanics parameter
kappa = zeros(n,1); % Fracture mechanics parameter

switch probtype
    case 'PlaneStress'
        
        for i = 1:n
            
            D{i} = [    1,  v(i),          0;   ...
                     v(i),     1,          0;   ...
                        0,     0, (1-v(i))/2    ].*(E(i)/(1-v(i)^2));
            
            kappa(i) = (3-v(i))/(1+v(i)); % in range [1.6667, 3]
            alpha(i) = 1;
            
        end
        
    case 'PlaneStrain'
        
        for i = 1:n
            
            D{i} = [ 1-v(i),    v(i),       0;  ...
                       v(i),  1-v(i),       0;  ...
                          0,       0, 0.5-v(i) ].*(E(i)/(1+v(i))/(1-2*v(i)));
            
            kappa(i) = 3-4*v(i); % in range [1.5, 3]
            alpha(i) = 1-v(i)^2; % in range [0.75, 1]
            
        end
        
    otherwise
        error('Problem type: ''PlaneStress'' or ''PlaneStrain'' ?')
        
end