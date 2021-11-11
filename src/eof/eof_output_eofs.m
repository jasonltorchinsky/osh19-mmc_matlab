function eof_output_eofs(params, eofs)

% Ensure the directory structure for output is created
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, 'eofs');

if ~(isfolder(out_path))
    mkdir(out_path);
end
if ~(isfolder(exp_path))
    mkdir(exp_path);
end
if ~(isfolder(component_path))
    mkdir(component_path);
end

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Get common parameters
[nx, ny, nz] = size(eofs.u_mjo1);
nz = nz - 1; % For consistency
nt = size(eofs.exp1, 2);

% Create output file, delete the current one if present
eofs_file_name = 'eofs.nc';
eofs_file = fullfile(component_path, eofs_file_name);
if isfile(eofs_file)
    delete(eofs_file)
end

Q_mode_code = 0;
if strcmpi(params.Q_mode, 'mid')
    Q_mode_code = 0;
elseif strcmpi(params.Q_mode, 'up')
    Q_mode_code = 1;
end

nccreate(eofs_file, 'Q_mode_code',...
    'Datatype', 'int64',...
    'Format', 'netcdf4');
ncwriteatt(eofs_file, 'Q_mode_code', 'description', ...
    'Mid- (0) or upper-mode (1) of moisture for EOF decomposition');
ncwriteatt(eofs_file, 'Q_mode_code', 'units', 'N/A');
ncwrite(eofs_file, 'Q_mode_code', Q_mode_code);

% Variances of EOFs used to get back to physical units
nccreate(eofs_file, 'u_std',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(eofs_file, 'u_std', 'description', ...
    'Standard deviation of zonal wind speed part of EOF to convert to physical units.');
ncwriteatt(eofs_file, 'u_std', 'units', 'm s^(-1)');
ncwrite(eofs_file, 'u_std', eofs.u_std);

nccreate(eofs_file, 'q_std',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(eofs_file, 'q_std', 'description', ...
    'Standard deviation of moisture part of EOF to convert to physical units.');
ncwriteatt(eofs_file, 'q_std', 'units', 'kg^1 kg^(-1)');
ncwrite(eofs_file, 'q_std', eofs.q_std);


% "Raw" EOFs, not scaled by variances
nccreate(eofs_file, 'raw_eof1',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'2x', 2*nx});
ncwriteatt(eofs_file, 'raw_eof1', 'description', ...
    '"Raw" first EOF, unscaled by variances');
ncwriteatt(eofs_file, 'raw_eof1', 'units', 'N/A');
ncwrite(eofs_file, 'raw_eof1', eofs.raw_eof1);

nccreate(eofs_file, 'raw_eof2',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'2x', 2*nx});
ncwriteatt(eofs_file, 'raw_eof2', 'description', ...
    '"Raw" second EOF, unscaled by variances');
ncwriteatt(eofs_file, 'raw_eof2', 'units', 'N/A');
ncwrite(eofs_file, 'raw_eof2', eofs.raw_eof2);

% "Raw" expansion coefficients
% Expansion coefficients
nccreate(eofs_file, 'raw_exp1',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'t', nt});
ncwriteatt(eofs_file, 'raw_exp1', 'description', ...
    'Unscaled expansion coefficients for first EOF');
ncwriteatt(eofs_file, 'raw_exp1', 'units', 'N/A');
ncwrite(eofs_file, 'raw_exp1', eofs.raw_exp1);

nccreate(eofs_file, 'raw_exp2',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'t', nt});
ncwriteatt(eofs_file, 'raw_exp2', 'description', ...
    'Unscaled expansion coefficients for second EOF');
ncwriteatt(eofs_file, 'raw_exp2', 'units', 'N/A');
ncwrite(eofs_file, 'raw_exp2', eofs.raw_exp2);


% EOFs
nccreate(eofs_file, 'u_eof1',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'x', nx});
ncwriteatt(eofs_file, 'u_eof1', 'description', ...
    'First EOF of zonal wind');
ncwriteatt(eofs_file, 'u_eof1', 'units', 'm s^(-1)');
ncwrite(eofs_file, 'u_eof1', eofs.u_eof1);

nccreate(eofs_file, 'u_eof2',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'x', nx});
ncwriteatt(eofs_file, 'u_eof2', 'description', ...
    'Second EOF of zonal wind');
ncwriteatt(eofs_file, 'u_eof2', 'units', 'm s^(-1)');
ncwrite(eofs_file, 'u_eof2', eofs.u_eof2);

nccreate(eofs_file, 'q_eof1',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'x', nx});
ncwriteatt(eofs_file, 'q_eof1', 'description', ...
    'First EOF of moisture');
ncwriteatt(eofs_file, 'q_eof1', 'units', 'kg kg^(-1)');
ncwrite(eofs_file, 'q_eof1', eofs.q_eof1);

nccreate(eofs_file, 'q_eof2',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'x', nx});
ncwriteatt(eofs_file, 'q_eof2', 'description', ...
    'Second EOF of moisture');
ncwriteatt(eofs_file, 'q_eof2', 'units', 'kg kg^(-1)');
ncwrite(eofs_file, 'q_eof2', eofs.q_eof2);

% Expansion coefficients
nccreate(eofs_file, 'exp1',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'t', nt});
ncwriteatt(eofs_file, 'exp1', 'description', ...
    'Expansion coefficients for first EOF');
ncwriteatt(eofs_file, 'exp1', 'units', 'N/A');
ncwrite(eofs_file, 'exp1', eofs.exp1);

nccreate(eofs_file, 'exp2',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'t', nt});
ncwriteatt(eofs_file, 'exp2', 'description', ...
    'Expansion coefficients for second EOF');
ncwriteatt(eofs_file, 'exp2', 'units', 'N/A');
ncwrite(eofs_file, 'exp2', eofs.exp2);

% Scaled covariance fraction
nccreate(eofs_file, 'scf',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'neof', 2*nx});
ncwriteatt(eofs_file, 'scf', 'description', ...
    'Scaled covariance fractions for EOFs');
ncwriteatt(eofs_file, 'scf', 'units', 'N/A');
ncwrite(eofs_file, 'scf', eofs.scf);

% Reconstructed MJOs
nccreate(eofs_file, 'u_mjo1',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz+1});
ncwriteatt(eofs_file, 'u_mjo1', 'description', ...
    'Zonal wind of reconstructed first MJO mode');
ncwriteatt(eofs_file, 'u_mjo1', 'units', 'm s^(-1)');
ncwrite(eofs_file, 'u_mjo1', eofs.u_mjo1);

nccreate(eofs_file, 'u_mjo2',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz+1});
ncwriteatt(eofs_file, 'u_mjo2', 'description', ...
    'Zonal wind of reconstructed second MJO mode');
ncwriteatt(eofs_file, 'u_mjo2', 'units', 'm s^(-1)');
ncwrite(eofs_file, 'u_mjo2', eofs.u_mjo2);

nccreate(eofs_file, 'q_mjo1',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz+1});
ncwriteatt(eofs_file, 'q_mjo1', 'description', ...
    'Moisture of reconstructed first MJO mode');
ncwriteatt(eofs_file, 'q_mjo1', 'units', 'kg kg^(-1)');
ncwrite(eofs_file, 'q_mjo1', eofs.q_mjo1);

nccreate(eofs_file, 'q_mjo2',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz+1});
ncwriteatt(eofs_file, 'q_mjo2', 'description', ...
    'Moisture of reconstructed second MJO mode');
ncwriteatt(eofs_file, 'q_mjo2', 'units', 'kg kg^(-1)');
ncwrite(eofs_file, 'q_mjo2', eofs.q_mjo2);
end

