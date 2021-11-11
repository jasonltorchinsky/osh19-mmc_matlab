function q = ens_unproj_q(dcm_params, dcm_grid, q_proj, eofs)

% Unpack some common parameters
nx = dcm_params.nx;
ny = dcm_params.ny;
nz = dcm_params.nz;

H = dcm_params.H;

yy = dcm_grid.yy;
zzW = dcm_grid.zzU;

dy = dcm_grid.dy;

% Scale meridional, vertical coordinate for basis functions
L = 1490; % Equatorial meridional length scale (km)
yy_norm = yy / L;
dy_norm = dy / L;
zzW_norm = pi * zzW / H;

% Get first baroclinic mode, zeroth parabolic cylinder function, standard
% deviation scale of u for projections
parab_cyl_0   = parab_cyl(yy_norm, 0);
q_clin_mode_1 = q_clin_mode(zzW_norm, 1);
q_clin_mode_2 = q_clin_mode(zzW_norm, 2);

Q_mode = eofs.Q_mode;
q_std = eofs.q_std;

q = zeros([ny, nx, nz + 1]);

% Redimensionalize q
q_proj_dim = q_proj * q_std;

% Unproject meridionally
q_mer = zeros([ny, nx]);
for ii = 1:nx
    q_mer(:,ii) = q_proj_dim(ii).*parab_cyl_0;
end

% Unproject vertically - NOTE: the unprojection is actually just the mid- or
% upper-barolicinic mode of moisture
if strcmpi(Q_mode, 'mid')
    for jj = 1:ny
        for ii = 1:nx
           q(jj, ii, :) = q_mer(jj, ii) * (q_clin_mode_1 - 1/4 * q_clin_mode_2); 
        end
    end
else
    for jj = 1:ny
        for ii = 1:nx
           q(jj, ii, :) = q_mer(jj, ii) * q_clin_mode_1; 
        end
    end
end

end

