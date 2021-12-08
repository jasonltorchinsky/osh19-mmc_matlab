function eofs_out = eof_recon_mjo(params, eofs_in)

eofs_out = eofs_in;

% Get path to output data.
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Read in run parameters
params_file = fullfile(component_path, 'params.nc');

nx = ncread(params_file, 'nx');
ny = ncread(params_file, 'ny');
nz = ncread(params_file, 'nz');

H  = ncread(params_file, 'H');

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

yy  = ncread(grid_file, 'yy');
zzU = ncread(grid_file, 'zzU');
zzW = ncread(grid_file, 'zzW');

dz  = ncread(grid_file, 'dz');

L = 1490; % Equatorial meridional length scale (km)
yy_norm  = yy / L;
zzU_norm = pi * (zzU + dz/2) / H;
zzW_norm = pi * zzW / H;

% Reconstruct the zonal wind part of the MJO
parab_cyl_0   = parab_cyl(yy_norm, 0);
u_clin_mode_1 = u_clin_mode(zzU_norm, 1);

eofs_out.u_mjo1 = zeros([ny, nx, nz + 1]);
eofs_out.u_mjo2 = zeros([ny, nx, nz + 1]);
for jj = 1:ny
    for ii = 1:nx
        for kk = 1:nz+1
            eofs_out.u_mjo1(jj, ii, kk) = eofs_in.u_eof1(ii) ...
                * parab_cyl_0(jj) * u_clin_mode_1(kk);
            eofs_out.u_mjo2(jj, ii, kk) = eofs_in.u_eof2(ii) ...
                * parab_cyl_0(jj) * u_clin_mode_1(kk);
        end
    end
end

% Reconstruct the moisture part of the MJO
Q_1 = zeros([ny, nx]);
Q_2 = zeros([ny, nx]);
for jj = 1:ny
    Q_1(jj, :) = eofs_in.q_eof1 .* parab_cyl_0(jj);
    Q_2(jj, :) = eofs_in.q_eof2 .* parab_cyl_0(jj);
end

% We get the q_1, q_2 of each moisture EOF by projecting Q onto the moisture
% baroclinic modes
q_clin_mode_1 = q_clin_mode(zzW_norm, 1);
q_clin_mode_2 = q_clin_mode(zzW_norm, 2);

Q_mode = params.Q_mode;

eofs_out.q_mjo1 = zeros([ny, nx, nz + 1]);
eofs_out.q_mjo2 = zeros([ny, nx, nz + 1]);
if strcmpi(Q_mode, 'mid')
    for kk = 1:nz+1
        eofs_out.q_mjo1(:,:,kk) = 1 / sqrt(2) * Q_1 .* (q_clin_mode_1(kk) + q_clin_mode_2(kk));
        eofs_out.q_mjo2(:,:,kk) = 1 / sqrt(2) * Q_2 .* (q_clin_mode_1(kk) + q_clin_mode_2(kk));
    end
else
    for kk = 1:nz+1
        eofs_out.q_mjo1(:,:,kk) = 1 / sqrt(2) * Q_1 .* (q_clin_mode_1(kk) - q_clin_mode_2(kk));
        eofs_out.q_mjo2(:,:,kk) = 1 / sqrt(2) * Q_2 .* (q_clin_mode_1(kk) - q_clin_mode_2(kk));
    end
end

end

