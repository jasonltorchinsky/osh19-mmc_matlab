function u = ens_unproj_u(dcm_params, dcm_grid, u_proj, eofs)

% Unpack some common parameters
nx = dcm_params.nx;
ny = dcm_params.ny;
nz = dcm_params.nz;

H = dcm_params.H;

yy = dcm_grid.yy;
zzU = dcm_grid.zzU;

dy = dcm_grid.dy;

% Scale meridional, vertical coordinate for basis functions
L = 1490; % Equatorial meridional length scale (km)
yy_norm = yy / L;
zzU_norm = pi * zzU / H;

% Get first baroclinic mode, zeroth parabolic cylinder function, standard
% deviation scale of u for projections
parab_cyl_0   = parab_cyl(yy_norm, 0);
u_clin_mode_1 = u_clin_mode(zzU_norm, 1);
u_std = eofs.u_std;

u = zeros([ny, nx, nz + 1]);

% Re-dimensionalize u_proj
u_proj_dim = u_proj * u_std;

% Unproject meridionally
u_mer = zeros([ny, nx]);
for ii = 1:nx
    u_mer(:,ii) = u_proj_dim(ii).*parab_cyl_0;
end

% Unproject vertically
for jj = 1:ny
    for ii = 1:nx
        u(jj, ii, :) = u_mer(jj, ii).*u_clin_mode_1;
    end
end

end

