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

ens_eofs.u_std = ncread(eof_file, 'u_std');
ens_eofs.q_std = ncread(eof_file, 'q_std');

ens_eofs.eof1 = ncread(eof_file, 'eof1');
ens_eofs.eof2 = ncread(eof_file, 'eof2');

ens_eofs.u_mjo1 = ncread(eof_file, 'u_mjo1');
ens_eofs.u_mjo2 = ncread(eof_file, 'u_mjo2');

ens_eofs.q_mjo1 = ncread(eof_file, 'q_mjo1');
ens_eofs.q_mjo2 = ncread(eof_file, 'q_mjo2');

Q_mode_code = ncread(eof_file, 'Q_mode_code');
if Q_mode_code == 0
    ens_eofs.Q_mode = 'mid';
elseif Q_mode_code == 1
    ens_eofs.Q_mode = 'up';
end
    

end

