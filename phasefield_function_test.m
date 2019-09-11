% material_phasefield_test.m
%  (273.15 + 23) * ((101325 + 2 * 1000 * 9.81) / 101325) ^ 0.2857 - 273.15
%
% Solve fracture problem in wood considering wood as a composite material
% made up of three material phases: earlywood (EW), transition wood, and 
% latewood (LW).
%
% Wood pine critical fracture energy:
%   Gc = 100 J/m^2
%
% Toughness of green spruce [Attack et al. 1961]: 
%   100 J/m^2 in the TR– direction,
%   180 J/m^2 in the RT direction
%
% Stress intensity factors for air-dried Douglas fir [Schniewind and Centeno (1973)] 
%   0.35 MPa.m^(1/2) for both RT and TR directions (i.e. no difference)
% 
% Assumptions
% -----------
%
% Early-wood mechanical properties:
% E = 40 MPa
% nu = 0.33
%
% Late-wood mechanical properties:
% E = 1200 MPa 
% nu = 0.33
%
% Critical fracture energy
% Gc = 100 J/m^23
% 
% =========================================================================

SIGMA = 20;
% SIGMA = 0;

L = 1;
H = 1;

xlim = [-L, 0];
ylim = [-H/2, H/2];

color_image = imread('material_phasefield.png');
grayscale_image = rgb2gray(color_image);

% Apply Gaussian filter
if SIGMA > 0
    grayscale_image = imgaussfilt(grayscale_image, SIGMA);
end

grayscale_image = im2double(grayscale_image);

[material_phasefield, ddx_material_phasefield, ddy_material_phasefield] = ...
    phasefield_function(grayscale_image, xlim, ylim);

x = linspace(xlim(1), xlim(2), 250);
y = linspace(ylim(1), ylim(2), 250);

[X, Y] = meshgrid(x,y);

X = X(:);
Y = Y(:);

fig_1 = figure(1);
imshow(grayscale_image)
axis equal

fig_2 = figure(2);
C = material_phasefield(X, Y);
scatter(X, Y, 3, C)
title('Material phasefield')
axis equal

fig_3 = figure(3);
C = ddx_material_phasefield(X, Y);
scatter(X, Y, 3, C)
title('d/dx')
axis equal

fig_4 = figure(4);
C = ddy_material_phasefield(X, Y);
scatter(X, Y, 3, C)
title('d/dy')
axis equal
