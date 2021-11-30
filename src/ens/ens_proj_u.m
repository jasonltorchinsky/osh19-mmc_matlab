function u_proj = ens_proj_u(dcm_params, dcm_grid, dcm_state, eofs)

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
dy_norm = dy / L;
zzU_norm = pi * zzU / H;

% Get first baroclinic mode, zeroth parabolic cylinder function, standard
% deviation scale of u for projections
parab_cyl_0   = parab_cyl(yy_norm, 0);
u_clin_mode_1 = u_clin_mode(zzU_norm, 1);
u_std = eofs.u_std;

% Get the normalization constants due to discretized norms
merid_norm = dy_norm * (parab_cyl_0 * parab_cyl_0.');
vert_norm  = (1/(nz+1)) * (u_clin_mode_1 * u_clin_mode_1.');

u_proj = zeros([1, nx]);

% Project u onto first baroclinic mode
u1 = zeros([ny, nx]);
for jj = 1:ny
    for ii = 1:nx
        u1(jj, ii) = (1/(nz+1)) ...
            * squeeze(u_clin_mode_1*squeeze(dcm_state.u(jj, ii, :))) ...
            * (1/vert_norm);
    end
end

% Project u1 onto zeroth parabolic cylinder function
for ii = 1:nx
    u_proj(1, ii) = dy_norm * squeeze(parab_cyl_0*u1(:, ii)) ...
        * (1/merid_norm);
end

% Scale by standard deviation to non-dimensionalize
u_proj = u_proj / u_std;

end

