% The code has been tested under Matlab R2015a. If older or newer versions
% are used there may be incompatibilities. Check them and adapt the code if
% required.

clearvars
clc
%% COMPONENTS CELL ARRAYS
% NOTE: Be sure to name the component identically to the .cosmo file name
alcohol = {'ethanol0.cosmo'};
alkane  = {'hexane.cosmo'};
ester   = {'methyloctanoate.cosmo'};

nalcohol = size(alcohol,1);
nalcano  = size(alkane, 1);
nester   = size(ester,  1);

%% CREATE A COSMO-RS JOB FOR A TERNARY SYSTEM CALCULATION
% (Most steps are shared for any other calculation)

% This creates a job structure that will be handled to the routines that
% compute the desired properties
job = new_cosmors();

% Make directory for new job (Prepare your own. Maybe:
% C/Users/user/Dekstop, or wherever you want to locate your working
% directory. In this directory is where the code will place .inp, .out, and
% .tab files. NOTE: do not follow the line below exactly. If you whish you
% can hard code your path
job.wdir = [ sprintf(job.wdir, 'J0029') 'HE/'];

% If the working directory does not exist, create it and add it to Matlab
% path
if ~isdir(job.wdir)    
    mkdir(job.wdir)    
end

addpath(job.wdir)

% This is the path to the .cosmo files database. This database must have
% the cosmotherm standard. This is, organized alphabetically into folders
% named by components' initial character (0,1,2...,a,b,...z). You may use
% built-in cosmotherm TZVP database
job.dbpath = 'D:/TERMO/TL DATA/DATABASE-COSMO/01_Organicos/02_BP-TZVP-COSMO';

% This is a header with calculation general settings
job.header2 = 'wtln unit=si notempty ehfile';

% Number of components in your mixtures
job.ncomponents  = 3;

% What kind of calculation you want to perform. In this case I selected
% <multinary> because there are more than two components in the mixture.
% Check cosmotherm user guide. In case you are working with binaries, just
% head the line with <binary>.
% Please, note that the line's xend and tk parameters are to be overwritten
% for each temperature-composition vector that is to be evaluated. This is
% why the text format (e.g. "%3.1f") is used instead of hard coding a
% certain value of composition and temperature.
job.jinstruction = 'multinary tk=%3.1f xstart={0 0 1} xend={%6.5f %6.5f %6.5f} xstep=2 HE_SPLIT';

% Number of .bat files in which to distribute the calculations
% This comes in handy when you're running many calculations because .bat
% files may be run in parallel. Inside each .bat file all calculations are
% run sequentially. For not too many calculations using just one bat file
% is enough though.
job.num_bat_files = 1;

% Create calculation files (.inp)
set_cosmors_inp_file(job, alcohol, alkane, ester );

%% COMPUTE EXCESS ENTHALPY FOR A TERNARY SYSTEM

% Creates a mesh of temperature and composition to evaluate the desired
% property
ptx = multimeshgrid(298,linspace(1e-9,1-1e-9,10),linspace(1e-9,1-1e-9,10));

% Remove rows with x1+x2>1
ptx(sum(ptx(:,2:3),2)>1,:)=[];

% Computes x3 as 1-x1-x2 and concatenates it to ptx matrix
ptx = [ptx 1-sum(ptx(:,2:3),2)];

% Get system's component names without extension
[~,compnames] = cellfun(@fileparts,...
                        [alcohol alkane ester],...
                        'uni',0);

% Define the current inp filename
job.system = [compnames{1} '+' compnames{2} '+' compnames{3} '-' job.calc];

% Compute excess enthalpy at each mesh point
[hxs_msh, mf_msh, hb_msh, vdw_msh] = hxs_cosmors( ptx, job );

%% HOW TO FINE TUNE A PARAMETER

% ---------------------------------------------------------------------
% Requirement 1: Modify job.header2 field including which COSMO-RS
% parameters are going to be optimized. (e.g., if we want to fine tune
% misfit energy contribution prefactor, cmfset, we shall do...).
% NOTE: 'wtln unit=si notempty ehfile' is common. What it is specific to
% this calculation is the 'cmfset=%3.4e'
job.header2 = 'wtln unit=si notempty ehfile cmfset=%3.4e';

% Requirement 2: Build a wrapper (here, optimization_wrapper_example.mat)
% function that will be called by the optimizer at every iteration. This
% wrapper function should modify job.header2 field with the new 
% parameter/s being tested.
%----------------------------------------------------------------------

% Optimizer call (here I use fminsearch, but any other optimization method
% could be used). 
% NOTE: This is a dummy example to demonstrate how fine tuning works. The
% variable "hxs_msh" has been calculated in the section above using
% cmfset=1. What I'm going to do is to set this parameter as 0.5 as an
% initial guess and see what results from the optimization process. Since
% evaluating the whole ptx matrix is computationally expensive, I'll use 
% just some rows (10) that I'll pick at random (see below).
x0 = 0.5;                          % Initial guess
idx = randsample(size(ptx,1), 10); % A sample of the mesh defined above to perform fittin demonstration
XDATA = ptx(idx,:);                % [T,x1,x2,x3]-experimental data.
YDATA = hxs_msh(idx,:);            % hE experimental data at XDATA.
[x,fval] = fminsearch( @(x) optimization_wrapper_example(x, XDATA, YDATA), x0,...
                       optimset('MaxIter',10) );
                   
% As expected, the final x (cmfset parameter) is approximately 1. So,
% everything is working great!