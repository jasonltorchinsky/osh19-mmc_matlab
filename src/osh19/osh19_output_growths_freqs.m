function osh19_output_growths_freqs(params, growths_freqs)

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
n_wavenums = growths_freqs.n_wavenums;
n_modes    = growths_freqs.n_modes;

% Create output file, delete the current one if present
growths_freqs_file_name = 'growths_freqs.nc';
growths_freqs_file = fullfile(component_path, growths_freqs_file_name);
if isfile(growths_freqs_file)
    delete(growths_freqs_file)
end

nccreate(growths_freqs_file, 'wavenums', ...
    'Datatype', 'double', ...
    'Format', 'netcdf4', ...
    'Dimensions', {'n_kx', n_wavenums});
ncwriteatt(growths_freqs_file, 'wavenums', 'description', ...
    'Zonal wavenumbers in initial condition');
ncwriteatt(growths_freqs_file, 'wavenums', 'units', 'N/A');
ncwrite(growths_freqs_file, 'wavenums', growths_freqs.wavenums);

nccreate(growths_freqs_file, 'modes', ...
    'Datatype', 'double', ...
    'Format', 'netcdf4', ...
    'Dimensions', {'n_modes', n_modes});
ncwriteatt(growths_freqs_file, 'modes', 'description', ...
    'Modes of zonal wavenumbers in initial condition');
ncwriteatt(growths_freqs_file, 'modes', 'units', 'N/A');
ncwrite(growths_freqs_file, 'modes', growths_freqs.modes);


nccreate(growths_freqs_file, 'growth_rates', ...
    'Datatype', 'double', ...
    'Format', 'netcdf4', ...
    'Dimensions', {'n_kx', n_wavenums, 'n_modes', n_modes});
ncwriteatt(growths_freqs_file, 'growth_rates', 'description', ...
    'Growth rates of each mode in each zonal wave number in initial condition');
ncwriteatt(growths_freqs_file, 'growth_rates', 'units', 'd^(-1)');
ncwrite(growths_freqs_file, 'growth_rates', growths_freqs.growth_rates);

nccreate(growths_freqs_file, 'freqs', ...
    'Datatype', 'double', ...
    'Format', 'netcdf4', ...
    'Dimensions', {'n_kx', n_wavenums, 'n_modes', n_modes});
ncwriteatt(growths_freqs_file, 'freqs', 'description', ...
    'Frequencies of each mode in each zonal wave number in initial condition');
ncwriteatt(growths_freqs_file, 'freqs', 'units', 'd^(-1)');
ncwrite(growths_freqs_file, 'freqs', growths_freqs.freqs);


end

