function rmse = optimization_wrapper_example(x, XDATA, YDATA, job)

% Update job.header2 instruction with the new-tested parameter value, x
job.header2 = sprintf(job.header2, x);

% Update system's inp file with the new header
thefile = fullfile(job.dir,job.system,'.inp');
% 1) Read old file into <A>
A = regexp( fileread( thefile ), '\n', 'split');

% 2) Replace file second line with updated job.header2
A{2} = job.header2;

% 3) Update old inp file to perform calculation
fid = fopen( thefile , 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid);

% Compute excess enthalpy at XDATA points with the updated COSMO-RS model
hxs = cosmors_hxs( XDATA, job );

% Compute some error metric between the data and the model. (Here, the root
% mean squared error is proposed)
rmse = sqrt( mean( (YDATA - hxs).^2 ) );