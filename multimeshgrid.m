% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is shared under the Creative Commons (CC BY-SA 4.0) license.
% Authorship of this piece of code is credited to Adriel Sosa PhD (email: 
% adriel.sosa@ulpgc.es) and Guillermo Alamo PhD
% The code is provided as is, so the author is not responsible for any kind 
% of flaw, inaccuracy or malfunctioning it may contain.
% Copyright (r) CC BY-SA 4.0, 2021 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IDX = multimeshgrid( varargin )
%MULTIMESHGRID Produces cartesian product matrix from n arrays
%   Useful to set all combinations of indices or elements of different
%   arrays. The output is not valid itself for 3D surfing, but the output
%   may be reshaped afterwards for this purpose.
%   Still, if the last were the case, its better to call meshgrid directly.


c = varargin;
L = cellfun(@(c) length(c),c);
LL = length(L);
nt = prod(L);
IDX = zeros(nt,LL);

for n=1:LL
    if size(c{n},1)==1
        c{n} = c{n}';
    end
    nr = prod(L(1:LL<n));
    np = prod(L(1:LL>n));
    IDX(:,n) = repmat(reshape(repmat(c{n}',[np 1]),[],1),[nr 1]);
end

                

