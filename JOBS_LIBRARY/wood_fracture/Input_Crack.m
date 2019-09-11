

%==========================================================================
% Crack Distribution
%==========================================================================


%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

nCrack = 0;  % number of cracks; nCrack = length(cCkCrd)
cCkCrd = []; % cell of crack coordinates: {[x1,y1;x2,y2;...];[];...}
mTpFrz = []; % crack tips to restrain (or freeze): [i_crk, i_tip;...]
mCkBox = []; % crack growth bounding box: [x_min,y_min;x_max,y_max],
             % i.e. confines cracks to smaller computational domain

%--------------------------------------------------------------------------
% Define cracks manually
%--------------------------------------------------------------------------

% % load cracks' file
% load([job_srcdir,'/Input_Crack/', ...
%     sprintf('job_crackID(%i).mat',job_crackID)]);

initial_crack_length = 1;

cCkCrd{1}  = [-initial_crack_length, 0; 1e-4, 0];

nCrack = length(cCkCrd);

%--------------------------------------------------------------------------
