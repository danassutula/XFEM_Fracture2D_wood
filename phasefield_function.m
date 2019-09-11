function [fh_F, fh_dFdx, fh_dFdy] = ...
    phasefield_function(grayscale_image, xlim, ylim)

if size(grayscale_image, 3) ~= 1
   error("`grayscale_image` must be a 2D matrix") 
end

I = im2double(grayscale_image);

nr = size(I, 1);
nc = size(I, 2);

x0 = xlim(1);
x1 = xlim(2);

y0 = ylim(1);
y1 = ylim(2);

L = x1 - x0;
H = y1 - y0;

dy = H / (nr - 1);
dx = L / (nc - 1); 

max_I = max(I, [], 'all');
min_I = min(I, [], 'all');

% Scale image intensity to range [0,1]
I(:) = (I - min_I) / (max_I - min_I);

% Compute image gradients
[dIdj, dIdi] = imgradientxy(I, 'sobel');

% NOTE: "prewitt" needs to be normalized by 6
% NOTE: "sobel" needs to be normalized by 8
% NOTE: "central" needs not be normalized

dIdx = dIdj / (8*dx);
dIdy =-dIdi / (8*dy);

function [i, j] = xy2sub(x, y)
   i = 1+round((y1-y(:))/dy);
   j = 1+round((x(:)-x0)/dx);
end

function value = F(x, y)
    [i, j] = xy2sub(x, y);
    ind = (j-1) * nr + i;
    value = I(ind);
end

function value = dFdx(x, y)
    [i, j] = xy2sub(x, y);
    ind = (j-1) * nr + i;
    value = dIdx(ind);
end

function value = dFdy(x, y)
    [i, j] = xy2sub(x, y);
    ind = (j-1) * nr + i;
    value = dIdy(ind);
end

fh_F = @F;
fh_dFdx = @dFdx;
fh_dFdy = @dFdy;

end
