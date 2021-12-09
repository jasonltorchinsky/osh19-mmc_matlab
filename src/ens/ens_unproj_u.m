function u = ens_unproj_u(dcm_params, dcm_grid, u_proj, eofs)

% Unpack some common parameters
nx = dcm_params.nx;
ny = dcm_params.ny;
nz = dcm_params.nz;

H = dcm_params.H;

yy = dcm_grid.yy;
zzU = dcm_grid.zzU;

dz = dcm_grid.dz;

% Scale meridional, vertical coordinate for basis functions
L = 1490; % Equatorial meridional length scale (km)
yy_norm = yy / L;
zzU_norm = pi * (zzU(2:nz+1) - dz/2) / (H - dz);

% Get first baroclinic mode, zeroth parabolic cylinder function, standard
% deviation scale of u for projections
parab_cyl_0   = parab_cyl(yy_norm, 0);
u_clin_mode_1 = u_clin_mode(zzU_norm, 1);
u_std = eofs.u_std/1000; % File is in m s^(-1), code is in km s^(-1)

u = zeros([ny, nx, nz + 1]);

% Dimensionalize, and unproject vertically and meridionally
for jj = 1:ny
    for ii = 1:nx
        for kk = 2:nz+1
            u(jj, ii, kk) = u_std * u_proj(ii) * parab_cyl_0(jj) ...
                * u_clin_mode_1(kk-1);
        end
    end
end

end

