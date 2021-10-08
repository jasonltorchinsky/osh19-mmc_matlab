function osh19_output_params(params)

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
params_file_name = 'params.nc';
params_file = fullfile(component_path, params_file_name);
if isfile(params_file)
    delete(params_file)
end

% Spatial, temporal grid spacing
nccreate(params_file, 'nx',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'nx', 'description', ...
    'Number of zonal grid points');
ncwriteatt(params_file, 'nx', 'units', 'N/A');
ncwrite(params_file, 'nx', params.nx);

nccreate(params_file, 'dy',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'dy', 'description', ...
    'Meridional grid spacing');
ncwriteatt(params_file, 'dy', 'units', 'km');
ncwrite(params_file, 'dy', grid.dy);

nccreate(params_file, 'dz',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'dz', 'description', ...
    'Vertical grid spacing');
ncwriteatt(params_file, 'dz', 'units', 'km');
ncwrite(params_file, 'dz', grid.dz);

nccreate(params_file, 'dt',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'dt', 'description', ...
    'Times-tep size');
ncwriteatt(params_file, 'dt', 'units', 's');
ncwrite(params_file, 'dt', grid.dt);

% Grids
nccreate(params_file, 'xx',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'x', nx});
ncwriteatt(params_file, 'xx', 'description', ...
    'Zonal grid points');
ncwriteatt(params_file, 'xx', 'units', 'km');
ncwrite(params_file, 'xx', grid.xx);

nccreate(params_file, 'yy',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny});
ncwriteatt(params_file, 'yy', 'description', ...
    'Meridional grid points');
ncwriteatt(params_file, 'yy', 'units', 'km');
ncwrite(params_file, 'yy', grid.yy);

nccreate(params_file, 'zzU',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'z', nz+1});
ncwriteatt(params_file, 'zzU', 'description', ...
    'Vertical grid levels for u, v, p');
ncwriteatt(params_file, 'zzU', 'units', 'km');
ncwrite(params_file, 'zzU', grid.zzU);

nccreate(params_file, 'zzW',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'z', nz+1});
ncwriteatt(params_file, 'zzW', 'description', ...
    'Vertical grid levels for w, theta, q');
ncwriteatt(params_file, 'zzW', 'units', 'km');
ncwrite(params_file, 'zzW', grid.zzW);

end

