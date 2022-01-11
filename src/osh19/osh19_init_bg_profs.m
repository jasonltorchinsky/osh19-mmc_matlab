function bg_profs = osh19_init_bg_profs(params, grid)

bg_profs = struct();

% Unpack some commonly used parameters
nx = params.nx;
ny = params.ny;
nz = params.nz;
H  = params.H;

dz = grid.dz;

% vertical grid for background moisture profile
zz = transpose(dz:dz:H-dz);
zz_full = transpose(0:dz:H);

% Unscaled baroclinic modes
sinz     = transpose(sin(zz * pi / H));
re_sinz  = reshape(sinz, [1, 1, length(sinz)]);
sinz_mat = repmat(re_sinz, [ny, nx, 1]);

sin2z     = transpose(sin(2 * zz * pi / H));
re_sin2z  = reshape(sin2z, [1, 1, length(sin2z)]);
sin2z_mat = repmat(re_sin2z, [ny, nx, 1]);


% Background moisture profiles
qbc_bot      = -params.B_vs * H;
qbc_top      = 0;
atemp        = qbc_bot / (1 - exp(-H / params.s_tilde));
btemp        = -atemp * exp(-H / params.s_tilde);
q_bg_vec_int = atemp * exp(-zz / params.s_tilde) + btemp;
q_bg_vec     = [qbc_bot; q_bg_vec_int; qbc_top];

ddz_q_bg_vec   = 1 / (2*dz) * (q_bg_vec(3:end) - q_bg_vec(1:end-2));
q_bg_merid_vec = ones(size(grid.yy)) ...
    - params.a * (1 - exp(-(grid.yy / params.L_tilde).^2/2));

[~, ~, q_bg_mat] = meshgrid(grid.xx, grid.yy, q_bg_vec);

for x_idx = 1:nx
    for z_idx = 1:nz+1
        q_bg_mat(:,x_idx,z_idx) = ...
            squeeze(q_bg_mat(:,x_idx,z_idx)).*transpose(q_bg_merid_vec);  % CHANGE 7/18/17!!
    end
end


% Background moisture damping profiles
tau_vec = params.tau_up + (params.tau_up - params.tau_mid) ...
    * (zz - H) / H;

D_h_vec = params.D_hUp + (params.D_hUp - params.D_hMid) ...
    * (zz - H) / H;

D_v_vec = params.D_v + (params.D_v - params.D_v) ...
    * (zz - H) / H;

% Background potential temperature profile
theta_bg_vec = params.B * zz_full;
[~, ~, theta_bg_mat] = meshgrid(grid.xx, grid.yy, theta_bg_vec);

% Save background profiles to bg_profs
bg_profs.zz       = zz;

bg_profs.sinz     = sinz;
bg_profs.sinz_mat = sinz_mat;

bg_profs.sin2z     = sin2z;
bg_profs.sin2z_mat = sin2z_mat;

bg_profs.q_bg_vec_int   = q_bg_vec_int;
bg_profs.q_bg_vec       = q_bg_vec;
bg_profs.q_bg_mat       = q_bg_mat;
bg_profs.ddz_q_bg_vec   = ddz_q_bg_vec;
bg_profs.q_bg_merid_vec = q_bg_merid_vec;

bg_profs.tau_vec = tau_vec;
bg_profs.D_h_vec = D_h_vec;
bg_profs.D_v_vec = D_v_vec;

bg_profs.theta_bg_mat = theta_bg_mat;

end

