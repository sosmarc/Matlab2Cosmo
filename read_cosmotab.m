% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is shared under the Creative Commons (CC BY-SA 4.0) license.
% Authorship of this piece of code is credited to Adriel Sosa PhD (email: 
% adriel.sosa@ulpgc.es)
% The code is provided as is, so the author is not responsible for any kind 
% of flaw, inaccuracy or malfunctioning it may contain.
% Copyright (r) CC BY-SA 4.0, 2021 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = read_cosmotab( ifile, gridsz )
%READ_COSMOTAB Import tab results into cell array of data
%   Detailed explanation goes here

if nargin<2
    gridsz = 29;
end

F = regexp( fileread(ifile), '\n', 'split');
F = F';

ini = double(cellfun(@(c) length(c)>=11 && ~isnan(str2double(c(1:11))),F));
ini(ini==0) = -1;
irow = find(diff(ini)==2)+1;

datsetsz = length(irow);
out = cell(datsetsz,1);

for m=1:datsetsz
    
    % Parse table header
    pos = column_name_end( F{irow(m)-1} );
    
    % Get formatted rows
    tmp = cellfun(@(r) format_row(r,pos),F(irow(m):(irow(m)+gridsz-1)), 'uni', 0);
        
    % Concatenate rows in tmp
    tmp = cat(1,tmp{:});   
    
    out{m} = cell2mat(tmp);
    
end

end

% Header preprocessing: To automate the reading of .tab files from
% COSMOtherm this function takes advantage of the fact that property names
% in table header columns are aligned to the last digit of the property
% values in remaining rows. So, if there is a whitespace in non-header rows 
% at the position right below that header's end-of-name char, it means that
% the value of that property in the inspected row is missing. That is
% useful to formate the read table.
% Hence, the functions below will look for the pattern ([^ ][ ]+) in the 
% header row to identify the positions where each column ends. This
% position vector will serve to explore remaining table rows to look for
% whitespaces and tag the missing properties in that row.
function  out = format_row( str, pos )

nprops = length(pos);

formatted_row = zeros(0, nprops );

% Get not missing positions
missing_pos = row_missing_pos( str, pos );
non_missing_pos = setdiff(1:nprops, missing_pos);

% Split the row by whitespaces. Collapse adjacent splitters
out = strsplit( str, ' ');

% Cast to double each splitted character
out = cellfun(@str2double, out, 'uni', 0);

% Detect right casting operations
cast_succeded = cell2mat( cellfun(@(v) ~isnan(v), out, 'uni', 0 ) );

% Keep succeeded castings
out = out( cast_succeded );

% Fill the formated row vector
formatted_row(non_missing_pos) = [out{:}];

% Zip the formatted vector into a cell
out = {formatted_row};

end

function pos = column_name_end(str)

pos = regexp(str, '([^ ][ ]+)' );
    
end

function val = row_missing_pos(str, pos)

val = regexp(str(pos), ' ');

end