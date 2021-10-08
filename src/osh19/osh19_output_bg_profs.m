function osh19_output_bg_profs(params, bg_profs)

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

% Unpack some common parameters
nx = params.nx;
ny = params.ny;
nz = params.nz;

% Create output file, delete the current one if present
bg_profs_file_name = 'bg_profs.nc';
bg_profs_file = fullfile(component_path, bg_profs_file_name);
if isfile(bg_profs_file)
    delete(bg_profs_file)
end


nccreate(bg_profs_file, 'q_bg_eq',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'z', nz + 1});
ncwriteatt(bg_profs_file, 'q_bg_eq', 'description', ...
    'Background moisture profile at equator');
ncwriteatt(bg_profs_file, 'q_bg_eq', 'units', 'kg kg^(-1)');
ncwrite(bg_profs_file, 'q_bg_eq', bg_profs.q_bg_vec);

nccreate(bg_profs_file, 'q_bg',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(bg_profs_file, 'q_bg', 'description', ...
    'Background moisture profile');
ncwriteatt(bg_profs_file, 'q_bg', 'units', 'kg kg^(-1)');
ncwrite(bg_profs_file, 'q_bg', bg_profs.q_bg_mat);

nccreate(bg_profs_file, 'theta_bg',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(bg_profs_file, 'theta_bg', 'description', ...
    'Background potential temperature profile');
ncwriteatt(bg_profs_file, 'theta_bg', 'units', 'K');
ncwrite(bg_profs_file, 'theta_bg', bg_profs.theta_bg_mat);

end

