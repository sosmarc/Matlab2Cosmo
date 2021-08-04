% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is shared under the Creative Commons (CC BY-SA 4.0) license.
% Authorship of this piece of code is credited to Adriel Sosa PhD (email: 
% adriel.sosa@ulpgc.es)
% The code is provided as is, so the author is not responsible for any kind 
% of flaw, inaccuracy or malfunctioning it may contain.
% Copyright (r) CC BY-SA 4.0, 2021 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ hxs,hmisf,hhb,hvdw ] = hxs_cosmors( ptx, cosmo )
%HXS_COSMORS Compute excess molar enthalpy according to COSMO-RS
%   The calculation is made at the Pressure-Temperature-Composition (ptx)
%   grid provided.
%   If Requested, the function also returns the contribution of Misfit,
%   Van der Waals and Hydrogen Bonds molecular interactions.

% Set the temperature vector
set_cosmors_pt(ptx,cosmo);

% Run calculation (Compute gamma at Tgrid)
dos(fullfile(cosmo.wdir, cosmo.batch));

% Get the number of components in the mixture
nc = size(ptx,2)-1;

% Get results
dat = read_cosmotab(fullfile(cosmo.wdir,[cosmo.system '.tab']),2);

dat = cat(1,dat{:});

% HXS + contributions
hxs   = dat(2:2:end,nc+1)*1000;

if nargout>1
    try
        hmisf = dat(2:2:end,(end-nc+1))*1000;
        hhb   = dat(2:2:end,(end-nc+2))*1000;
        hvdw  = dat(2:2:end,(end-nc+3))*1000;
    catch
        hmisf = [];
        hhb   = [];
        hvdw  = [];
    end
end