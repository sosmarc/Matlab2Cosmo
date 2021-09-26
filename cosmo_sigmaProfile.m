% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is shared under the Creative Commons (CC BY-SA 4.0) license.
% Authorship of this piece of code is credited to Adriel Sosa PhD (email: 
% adriel.sosa@ulpgc.es)
% The code is provided as is, so the author is not responsible for any kind 
% of flaw, inaccuracy or malfunctioning it may contain.
% Copyright (r) CC BY-SA 4.0, 2021 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sprof, surfSigInfoMat, cosmoSigma] = cosmo_sigmaProfile( comp, job, sigmaBinWidth, sigmaBounds, av )

if nargin<5 || isempty(av)
    av = 0.78389;   % Averaging surface in A^2 - Default COSMO parameter   (According to Sandler, this would be: 2.1002755 A)
                    %                                                      (According to Klamt, this would be: 0.78389 A)     
end
if nargin<4 || isempty(sigmaBounds)
    sigmaBounds = [-3 3];    % e·nm^-2
end
if length(sigmaBounds)~=2 || diff(sigmaBounds)<0
    error('COSMO:sigmaProfArgs', 'Wrong sigmaBounds. Argument must have length 2 and upper bound must be greater than lower bound');
end

if nargin<3 || isempty(sigmaBinWidth)
    sigmaBinWidth = 0.1;
end
if nargin < 2
    error('COSMO:sigmaProfArgs', 'Not enough input arguments, <2');
end

% Load .cosmo file and format segment-surface info matrix
surfSigInfoMat = get_sigSurfInfoMat( comp, job );

% Average surface-element charge-density (see Klamt et al., 1998. J. Phys.
% Chem. A, 102, 5074-5085)
cosmoSigma = sigma_averaging( surfSigInfoMat, av );

% Extract segment surface (in A^2)
sigmaSurf  = surfSigInfoMat(:,7);

% Get sigma profile (Surface-charge density from e·A^-2 to e·nm^-2 )
sprof = getSigmaProfile( cosmoSigma * 100, sigmaSurf, sigmaBinWidth, sigmaBounds );

% Smooth sigma-profile in case of Klamt-averaging
isKlamt = abs(av/0.78389-1)<1e-5;

if isKlamt
    sprof = klamtPostProcess( sprof );
end

end

function row = sigmaProfFomatRow( row )

if isempty(row)
    return
end

% Trim leading and trailing whitspaces
row = strtrim( row );

% Try to split by whitespaces
row = strsplit(row, ' ');

% If the previous action does not work try splitting by tabs
if length(row)==1
    row = strsplit(row, '\t');
elseif length(row)~=9
    error( 'COSMO:sigmaProfParser',...
           'It seems that read file is not a .cosmo file or has wrong format' );
end

% Convert
try
row = cellfun( @str2double, row, 'uni', 0);
catch ex
    disp('Error parsing sigma surface info. Data could not be coerced to double.');
    rethrow(ex);
end

row = cell2mat(row);

end

function sigmaProf = getSigmaProfile( cosmoSigma, A, sigmaBW, sigmaBounds )

Nsigma        = floor(diff(sigmaBounds)/sigmaBW)+1;
sigmaProf     = zeros(Nsigma,2);
sigmaProf(:,1)= sigmaBounds(1):sigmaBW:sigmaBounds(2);
sigmaPosReal  = (cosmoSigma - sigmaBounds(1))/sigmaBW;
sigmaPosInt   = floor(sigmaPosReal)+1;

% Sigma profile is made as the cummulative surface of the whole
% molecule with charge density equal to sigma. Because of the
% discretization of the sigma-domain, each surface-element area is assigned
% proportionally to two consecutive grid points of the sigma profile. See
% below
for i=1:(Nsigma-1)
    idx = find(sigmaPosInt==i);
    sigmaProf(i,2)   = sigmaProf(i,2)  +sum( A(idx).*( sigmaProf(i+1,1)-cosmoSigma(idx) )/sigmaBW );
    sigmaProf(i+1,2) = sigmaProf(i+1,2)+sum( A(idx).*( cosmoSigma(idx)-sigmaProf(i,1) )/sigmaBW ); 
end

end

function sigmaProf = klamtPostProcess( sigmaProf )
% In case Klamt averaging is being applied, smooth the calculated sigma profile
% using a Savitzky-Golay filter of window size 3

zeroPaddedSigmaProf = [0; sigmaProf(:,2); 0];
for i=2:(size(sigmaProf,1)+1)
    sigmaProf(i-1,2) = sum(zeroPaddedSigmaProf([i-1 i i+1])/3);
end

end

function sigSurfInfoMat = get_sigSurfInfoMat( comp, job )

% Ensure that the component is in the db
cap = comp(1);
if ~any( regexp(cap,'[a-zA-Z1-9]','ONCE')==1 )
    cap = num2str(0);
end
comp_formatted_name = fullfile(job.dbpath, lower(cap), [comp '.cosmo'] );

if ~exist(comp_formatted_name,'file')
    error('COSMO:sigmaProfReading',...
          'Provided component is not present in the provided database');
end

% Load component's cosmo file
try
    A = regexp( fileread(comp_formatted_name), '\n', 'split');
    A = A';
catch ex
    error('COSMO:sigmaProfReading',...
          ['Provided component is not present in the provided database.'...
           'Subyacent error:\n', ex.message]);
end

% Find <segment_information> tag
init = find( cell2mat( cellfun(@(str) strcmpi( strtrim(str),'$segment_information'),...
                               A, 'uni',0 ) ) );                              
init = init + 11;    % Apply offset to start reading surface segment info

% Get surface segment info matrix
sigSurfInfoMat = cellfun( @(row) sigmaProfFomatRow(row), A(init:end),...
                          'uni',0 );


% Concatenate (numeric) rows
NumericRows = cellfun(@isnumeric, sigSurfInfoMat);
sigSurfInfoMat = cat(1,sigSurfInfoMat{NumericRows});

end

function s_corr = sigma_averaging( sigSurfInfoMat, av )

N    = size( sigSurfInfoMat, 1 );          % Number of surface-segments (SS)
r2av = av / pi;                            % Averaging squared-radius
s_mu  = sigSurfInfoMat(:,8);               % surface-segment charge-density
s_corr = zeros(size(sigSurfInfoMat,1),1);  % Corrected SSCD

% Convert surface-segment positions form a.u to Amstrongs
sigSurfInfoMat(:,3:5) = sigSurfInfoMat(:,3:5) / 1.8897259885789;

% Compute squared-segment radius
r2mu = sigSurfInfoMat(:,7) / pi;

% Compute corrected surface-segment charge-density (SSCD)
for v=1:N    
    pos_v  = sigSurfInfoMat(v,3:5);
    d2_v_mu = sum( ( sigSurfInfoMat(:,3:5) - ones(N,1)*pos_v ).^2, 2);
    
    s_corr(v) = sum(s_mu.* (r2mu.*r2av)./(r2mu+r2av) .* exp(- d2_v_mu./ (r2mu+r2av) ) ) / ...
                sum((r2mu.*r2av)./(r2mu+r2av) .* exp(- d2_v_mu./ (r2mu+r2av) ) );
end

end
