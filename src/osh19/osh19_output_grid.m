function osh19_output_grid(params, grid)

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
grid_file_name = 'grid.nc';
grid_file = fullfile(component_path, grid_file_name);
if isfile(grid_file)
    delete(grid_file)
end

% Spatial, temporal grid spacing
nccreate(grid_file, 'dx',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(grid_file, 'dx', 'description', ...
    'Zonal grid spacing');
ncwriteatt(grid_file, 'dx', 'units', 'km');
ncwrite(grid_file, 'dx', grid.dx);

nccreate(grid_file, 'dy',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(grid_file, 'dy', 'description', ...
    'Meridional grid spacing');
ncwriteatt(grid_file, 'dy', 'units', 'km');
ncwrite(grid_file, 'dy', grid.dy);

nccreate(grid_file, 'dz',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(grid_file, 'dz', 'description', ...
    'Vertical grid spacing');
ncwriteatt(grid_file, 'dz', 'units', 'km');
ncwrite(grid_file, 'dz', grid.dz);

nccreate(grid_file, 'dt',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(grid_file, 'dt', 'description', ...
    'Times-tep size');
ncwriteatt(grid_file, 'dt', 'units', 's');
ncwrite(grid_file, 'dt', grid.dt);

% Grids
nccreate(grid_file, 'xx',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'x', nx});
ncwriteatt(grid_file, 'xx', 'description', ...
    'Zonal grid points');
ncwriteatt(grid_file, 'xx', 'units', 'km');
ncwrite(grid_file, 'xx', grid.xx);

nccreate(grid_file, 'yy',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny});
ncwriteatt(grid_file, 'yy', 'description', ...
    'Meridional grid points');
ncwriteatt(grid_file, 'yy', 'units', 'km');
ncwrite(grid_file, 'yy', grid.yy);

nccreate(grid_file, 'zzU',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'z', nz+1});
ncwriteatt(grid_file, 'zzU', 'description', ...
    'Vertical grid levels for u, v, p');
ncwriteatt(grid_file, 'zzU', 'units', 'km');
ncwrite(grid_file, 'zzU', grid.zzU);

nccreate(grid_file, 'zzW',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'z', nz+1});
ncwriteatt(grid_file, 'zzW', 'description', ...
    'Vertical grid levels for w, theta, q');
ncwriteatt(grid_file, 'zzW', 'units', 'km');
ncwrite(grid_file, 'zzW', grid.zzW);

end

