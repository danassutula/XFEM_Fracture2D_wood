
%--------------------------------------------------------------------------
% Material mechanical parameters
%--------------------------------------------------------------------------

% 2D assumtions
problemType = 'PlaneStrain'; % 'PlaneStress' or 'PlaneStrain'

% Units (plotting)
lengthUnits = 'mm';

nPhase = 2;

E = zeros(nPhase,1);
v = zeros(nPhase,1);

% Late wood
E(1) = 1200; % MPa (or MN/m^2, or N/mm^2)
v(1) = 0.33;

% Early wood
E(2) = 40; % MPa (or MN/m^2, or N/mm^2)
v(2) = 0.33;

% Material toughness
K_crt = 1; % MPa.m^(1/2)

% Material fracture energy
G_crt = 100e-3; % J/m^2 == 1e-3 N/mm 

%--------------------------------------------------------------------------
% Fracture mechanics parameters
%--------------------------------------------------------------------------

E_nearcracktip = mean(E);
v_nearcracktip = mean(v);

G_nearcracktip = 0.5*E_nearcracktip/(1+v_nearcracktip);

switch problemType
    case 'PlaneStress'
        
        kappa = (3-v_nearcracktip)/(1+v_nearcracktip); % in range [1.6667, 3]
        alpha = 1;
            
    case 'PlaneStrain'

        % fracture parameters
        kappa = 3-4*v_nearcracktip; % in range [1.5, 3]
        alpha = 1-v_nearcracktip^2; % in range [0.75, 1]

    otherwise
        error('Problem type: ''PlaneStress'' or ''PlaneStrain'' ?')
        
end

E_nearcracktip_star = E_nearcracktip / alpha;

%--------------------------------------------------------------------------
% Interpolation
%--------------------------------------------------------------------------

SIGMA = 0;

domain_xlim = [-5,0];
domain_ylim = [-2.00,2.00];

E = E(:);
v = v(:);

[cDMatx, ~, ~] = MtxConstit(nPhase, E, v, problemType);

% FIXME: Need to circumvent an error check in `Auxiliary_Discrete`
nPhase = 1;

% The easiest way to modify the program in order to allow for a function-
% like constitutive matrix is to define the constitutive-matrix function 
% as a global variable.

global constitutive_field % Constitutive relation as function
global constitutive_field_ddn % Directional derivative as function

C0 = cDMatx{1};
C1 = cDMatx{2};
dC = C1 - C0;

material_phasefield_image = rgb2gray(imread('material_phasefield.png'));

if SIGMA > 0
    material_phasefield_image = imgaussfilt(material_phasefield_image, SIGMA);
end

[fh_F, fh_dFdx, fh_dFdy] = ...
    phasefield_function(material_phasefield_image, domain_xlim, domain_ylim);

[constitutive_field, constitutive_field_ddn] = ...
    constitutive_interpolation(C0, C1, fh_F, fh_dFdx, fh_dFdy);

% constitutive_field = @(x) C0 + dC .* (-x(:,1)/5);
% constitutive_field_ddn = @(x,n) dC .* (-1/5);

% x = linspace(domain_xlim(1), domain_xlim(2), 350);
% y = linspace(domain_ylim(1), domain_ylim(2), 350);
% 
% [X, Y] = meshgrid(x,y);
% 
% X = X(:);
% Y = Y(:);
% 
% % fig_1 = figure(1);
% % imshow(material_phasefield_image)
% % axis equal
% 
% C = zeros(length(X),1); 
% dCdx = zeros(length(X),1); 
% dCdy = zeros(length(X),1); 
% 
% for i = 1:length(C)
%     C(i) = fh_F(X(i), Y(i));
%     dCdx(i) = fh_dFdx(X(i), Y(i));
%     dCdy(i) = fh_dFdy(X(i), Y(i));
% end
% 
% fig_2 = figure(51);
% scatter(X, Y, 3, C)
% title({'Material phasefield',...
%        sprintf('[Gaussian filtered with sigma = %d]',SIGMA)})
% colorbar
% axis equal

