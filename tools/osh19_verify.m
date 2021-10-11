function osh19_verify(params)

% Get path to output data.
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Read in run parameters
params_file = fullfile(component_path, 'params.nc');

sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

out_idxs = 0:1:floor(sim_days / out_freq);

% Read in the data from the verified run
u_ver = load('utot_verified.mat').utot_verified;

% Read data in for each time-step from the verified run and the default run of
% the code.
for out_idx = out_idxs
    
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    % Get total zonal wind
    u = ncread(state_file, 'u');
    
    disp(max(abs(u_ver(:,:,:,out_idx+1) - u(:,:,:)/1000), [], 'all'));
    % Note: Verified data has wind-speeds in m/s.
    
end

end

