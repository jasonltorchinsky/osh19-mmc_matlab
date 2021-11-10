function ens_eofs = ens_read_eofs(ens_params)

ens_eofs = struct();
eof_params = ens_params.eof_params;

% Ensure the directory structure for output is present
out_path       = eof_params.out_path;
exp_path       = fullfile(out_path, eof_params.exp_name);
component_path = fullfile(exp_path, eof_params.component_name);
eof_path       = fullfile(exp_path, 'eofs');

addpath(out_path);
addpath(exp_path);
addpath(component_path);
addpath(eof_path);

% Read in raw EOFs, variances needed to re-dimensionalize them later.
eof_file_name = 'eofs.nc';
eof_file = fullfile(eof_path, eof_file_name);

ens_eofs.u_var = ncread(eof_file, 'u_var');
ens_eofs.q_var = ncread(eof_file, 'q_var');

ens_eofs.raw_eof1 = ncread(eof_file, 'raw_eof1');
ens_eofs.raw_eof2 = ncread(eof_file, 'raw_eof2');

end

