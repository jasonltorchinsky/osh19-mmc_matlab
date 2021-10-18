function cmg14_output_state(params, state, out_idx, time)

% Ensure the directory structure for output is created
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

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


% Create output file, delete the current one if present
state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
state_file = fullfile(component_path, state_file_name);
if isfile(state_file)
    delete(state_file)
end

% Current time in simulation
nccreate(state_file, 't',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(state_file, 't', 'description', ...
    'Time of current state');
ncwriteatt(state_file, 't', 'units', 's');
ncwrite(state_file, 't', time);

% State variables
nccreate(state_file, 'u_1',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(state_file, 'u_1', 'description', ...
    'MJO 1');
ncwriteatt(state_file, 'u_1', 'units', 'N/A');
ncwrite(state_file, 'u_1', state.u_1);

nccreate(state_file, 'u_2',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(state_file, 'u_2', 'description', ...
    'MJO 2');
ncwriteatt(state_file, 'u_2', 'units', 'N/A');
ncwrite(state_file, 'u_2', state.u_2);

nccreate(state_file, 'v',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(state_file, 'v', 'description', ...
    'Stochastic damping');
ncwriteatt(state_file, 'v', 'units', 'N/A');
ncwrite(state_file, 'v', state.v);

nccreate(state_file, 'w_u',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(state_file, 'w_u', 'description', ...
    'Stochastic phase');
ncwriteatt(state_file, 'w_u', 'units', 'N/A');
ncwrite(state_file, 'w_u', state.w_u);

end

