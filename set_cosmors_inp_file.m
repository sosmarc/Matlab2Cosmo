% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is shared under the Creative Commons (CC BY-SA 4.0) license.
% Authorship of this piece of code is credited to Adriel Sosa PhD (email: 
% adriel.sosa@ulpgc.es)
% The code is provided as is, so the author is not responsible for any kind 
% of flaw, inaccuracy or malfunctioning it may contain.
% Copyright (r) CC BY-SA 4.0, 2021 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status,varargout] = set_cosmors_inp_file(job,varargin)
% SET_COSMORS_INP_FILE Sets inp files with desired calculation
%    Generates as many inp files as combinations can be made with provided
%    pure component names in varargin. varargin is, in fact, a cellarray of
%    cellstrings.
%    Pure component names must match those of the .cosmo files since the
%    code seeks them by name.

% Check if varargin is not empty
if isempty(varargin)
    error('At least one cellstr with component names must be provided')
end

% Check if varargin are strcells
if ~all( cellfun(@iscellstr, varargin, 'uni', 1) )
    error('Inputs must be cellstr objects with names of components in the mixture')
else
    comp_sets = varargin;
end

num_comps_by_set = cellfun(@(set) sum( 1 - cell2mat(cellfun(@(x) any(isnan(x)),set,'uni',0) ) ),...
                           comp_sets, 'uni', 1);
nc = length( num_comps_by_set );

% Number of combinations
numcases = prod(num_comps_by_set);
seqs = arrayfun(@(n) linspace(1,n,n), num_comps_by_set, 'uni', 0);
idx = multimeshgrid( seqs{:} );

% COSMO script files---------------------------------------------------
% A) 1 .inp file per system
% B) Packing calculus in bat files
jobs_per_bat = repmat(floor(numcases/job.num_bat_files),[job.num_bat_files,1]);
modd = numcases - sum(jobs_per_bat);
jobs_per_bat(job.num_bat_files:-1:(job.num_bat_files-modd+1)) = ...
                        jobs_per_bat(job.num_bat_files:-1:(job.num_bat_files-modd+1)) + 1;


% Write bat files
status = 0;
try
    jcount = 0;
    for m=1:job.num_bat_files
        
        % Open bat file
        bfilename = fullfile(job.wdir,sprintf('run_batch_B%07.0f.bat',...
            m)...
            );
        fbat= fopen(bfilename,'w');
        
        % Tell cmd to remain open after exectution
        fprintf(fbat,'@echo off\n');
        fprintf(fbat,sprintf('cd "%s"\n',job.wdir));
        
        % Write input files belongig to m-esim batch file
        for n = 1:jobs_per_bat(m)
            
            % Increase counter value
            jcount = jcount + 1;
            
            % Component and cosmo-file names in current job
            system_compnames = cell(1,nc);
            cosmonames       = cell(1,nc);

            % Loop over the number of compounds in the mixture
            for i=1:nc
                cap = comp_sets{i}{idx(jcount,i)}(1);
                if ~any( regexp(cap,'[a-zA-Z1-9]','ONCE')==1 )
                    cap = num2str(0);
                end
                cosmonames{i} = sprintf(job.comp, comp_sets{i}{idx(jcount,i)}, job.dbpath, lower(cap));
                
                [~,tmp] = fileparts( comp_sets{i}{idx(jcount,i)} );
                system_compnames{i} = tmp;
            end
            
            % Create the job (.inp) filename
            inp_file_pat = [repmat('%s+',[1 nc-1]) '%s-%s.inp'];
            jfilename = fullfile(job.wdir,sprintf(inp_file_pat,...
                                                  system_compnames{:},...
                                                  job.calc) );
            
            % Open job filename to define calculation
            fjob = fopen(jfilename,'w');
            
            % Headers
            fprintf(fjob,'%s\n',job.header1);
            fprintf(fjob,'%s\n',job.header2);
            fprintf(fjob,'!!------ Datos de la mezcla -------!!\n');
                        
            % Components
            for i=1:nc
                fprintf(fjob,cosmonames{i});
            end
            
            % Instructions
            for l = 1:job.Tgridsz(jcount)
                fprintf(fjob,'%s\n',job.jinstruction);
            end
            
            % Close current job file
            fclose(fjob);
            
            % Drop it into current bat file
            [~,jn,jxt] = fileparts(jfilename);
%             fprintf(fbat,'echo Processing %d file from %d in .bat file...\n',...
%                 n, jobs_per_bat(m));
            fprintf(fbat,job.binstruction,[jn jxt]);
        end
        
        % Close cmd window after run
        fprintf(fbat,'exit');
        
        % Close current bat file
        fclose(fbat);
    end
catch
    status = -1;
end
