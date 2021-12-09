function q_proj = ens_proj_q(dcm_params, dcm_grid, q, eofs)

% Unpack some common parameters
nx = dcm_params.nx;
ny = dcm_params.ny;
nz = dcm_params.nz;

H = dcm_params.H;

yy = dcm_grid.yy;
zzW = dcm_grid.zzW;

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

% Get the normalization constants due to discretized norms
merid_norm  = dy_norm * (parab_cyl_0 * parab_cyl_0.');
vert_norm_1 = (1/(nz+1)) * (q_clin_mode_1 * q_clin_mode_1.');
vert_norm_2 = (1/(nz+1)) * (q_clin_mode_2 * q_clin_mode_2.');

q_proj = zeros([1, nx]);

% Project q onto first two baroclinic modes
q1 = zeros([ny, nx]);
q2 = zeros([ny, nx]);
for jj = 1:ny
    for ii = 1:nx
        q1(jj, ii) = (1/(nz+1)) ...
            * squeeze(q_clin_mode_1*squeeze(q(jj, ii, :))) ...
            * (1/vert_norm_1);
        q2(jj, ii) = (1/(nz+1)) ...
            * squeeze(q_clin_mode_2*squeeze(q(jj, ii, :))) ...
            * (1/vert_norm_2);
    end
end

% Pick Q_mid or Q_up REFURB, ACTUALLY LOW
if strcmpi(Q_mode, 'mid') % is actually low here
    Q = 1 / sqrt(3) * (q1 * q_clin_mode(pi/3, 1) - q2 * q_clin_mode(pi/3, 2));
else
    Q = 1 / sqrt(3) * (q1 * q_clin_mode(2*pi/3, 1) + q2 * q_clin_mode(2*pi/3, 2));
end

% % Pick Q_mid or Q_up
% if strcmpi(Q_mode, 'mid')
%     Q = q1 * sqrt(2) * sin(pi/2);
% else
%     Q = q1 * sqrt(2) * sin(2*pi/3) + q2 * 2 * sqrt(2) * sin(4*pi/3);
% end

% Project Q onto zeroth parabolic cylinder function
for ii = 1:nx
    q_proj(1, ii) = dy_norm * squeeze(parab_cyl_0*Q(:, ii)) ...
        * (1/merid_norm);
end

% Scale by standard deviation to non-dimensionalize
q_proj = q_proj / q_std;

end

