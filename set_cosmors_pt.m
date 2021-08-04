% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is shared under the Creative Commons (CC BY-SA 4.0) license.
% Authorship of this piece of code is credited to Adriel Sosa PhD (email: 
% adriel.sosa@ulpgc.es)
% The code is provided as is, so the author is not responsible for any kind 
% of flaw, inaccuracy or malfunctioning it may contain.
% Copyright (r) CC BY-SA 4.0, 2021 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function status = set_cosmors_pt( ptx, cosmo )
%SET_COSMORS_PT Modify Pressure-Temperature-Composition conditions in .inp
%   Detailed explanation goes here

status = 1;
try
    [uT,ncol] = size(ptx);
    ptx = cellfun(@(row) num2cell(row),...
                  mat2cell(ptx, ones(uT,1), ncol),...
                  'uni',0);
    
    fname = fullfile(cosmo.wdir,[cosmo.system '.inp']);
    % Get old file
    A = regexp( fileread(fname), '\n', 'split');
    rows2skip = 3 + ncol;
    A(rows2skip:(rows2skip-1+uT)) = cellfun(@(row) sprintf(cosmo.jinstruction,row{:}),...
                                            ptx,...
                                           'uni', 0);
    
    % Set new file
    fid = fopen(fname, 'w');
    fprintf(fid, '%s\n', A{:});
    fclose(fid);
catch ex
    disp(ex.message)
    status = 0;
end

