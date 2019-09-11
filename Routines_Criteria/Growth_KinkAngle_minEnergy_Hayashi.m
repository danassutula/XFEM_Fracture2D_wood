function b = Growth_KinkAngle_minEnergy_Hayashi(r)

% Arbitrary kink angle. Solution is interpolated using values at points.

% r = SIF ratio, K2/K1
% b = crack tip kink angle.

tol = 1e-4;

rng_b = -acos(1/3)-(-75.8/180*pi); % <6/180*pi; max difference for pure mode-II
dff_b = 0.1/180*pi;
abs_r = abs(r);

w = pi*[0.0 -0.04 -0.08 -0.12 -0.16 -0.20 -0.24 -0.28 -0.32 -0.36 -0.40 -0.44 -0.48 -0.52 -0.56 -0.60 -0.64 -0.68 -0.72 -0.76 -0.80];

kI1 = [1.0 0.99410 0.97655 0.94794 0.90913 0.86127 0.80579 0.74427 0.67837 0.60981 0.54024 0.47126 0.40426 0.34049 0.28095 0.22632 0.17664 0.12248 0.10537 0.07269 0.04716];
kI2 = [0.0 0.18770 0.37069 0.54441 0.70469 0.84784 0.97083 1.07134 1.14784 1.19960 1.22672 1.22996 1.21082 1.17129 1.11390 1.04149 0.95746 0.87185 0.75787 0.65584 0.55040];

kII1 = [0.0 -0.06251 -0.12320 -0.18029 -0.23219 -0.27751 -0.31514 -0.34430 -0.36449 -0.37560 -0.37782 -0.37159 -0.35768 -0.33705 -0.31080 -0.28024 -0.24697 -0.21803 -0.17069 -0.13673 -0.10410];
kII2 = [1.0  0.98772  0.95131  0.89211  0.81224  0.71460  0.60262  0.48016  0.35137  0.22048  0.09165 -0.03116 -0.14440 -0.24495 -0.33023 -0.39818 -0.44724 -0.47259 -0.49145 -0.48255 -0.45761];

ke = sqrt((kI1+kI2*abs_r).^2+(kII1+kII2*abs_r).^2);

pp = spline(w,ke); % cubic spline interpolation

b0 = Growth_KinkAngle_maxTension(abs_r); % always lesser

bi = linspace(b0-rng_b,b0,50); % db ~ 0.1 (deg.)
[~,i] = max(ppval(pp,bi));
b = sign(r)*bi(i);

end