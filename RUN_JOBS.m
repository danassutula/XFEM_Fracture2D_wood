
%==========================================================================
% RUN JOBS
%==========================================================================

clear all;
close all;

job_titles = { ...
    'wood_fracture',
};

for job_title_i = 1:length(job_titles)
    job_title = job_titles{job_title_i};
    
    RUN_JOB
    
end

% exit(0)
