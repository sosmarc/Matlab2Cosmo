% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is shared under the Creative Commons (CC BY-SA 4.0) license.
% Authorship of this piece of code is credited to Adriel Sosa PhD (email: 
% adriel.sosa@ulpgc.es)
% The code is provided as is, so the author is not responsible for any kind 
% of flaw, inaccuracy or malfunctioning it may contain.
% Copyright (r) CC BY-SA 4.0, 2021 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lgamma = gamma_cosmors( ptx, cosmo )
%GAMMA_COSMORS Compute log activity coefficients according to COSMO-RS
%   The calculation is made at the Pressure-Temperature-Composition (ptx)
%   grid provided.

% Set the temperature vector
set_cosmors_pt(ptx,cosmo)

% Run calculation (Compute gamma at Tgrid)
dos(fullfile(cosmo.wdir, cosmo.batch));

% Get results
dat = read_cosmotab(fullfile(cosmo.wdir,[cosmo.system '.tab']),2);

dat = cat(1,dat{:});
lgamma = dat(2:2:end,8:9);