function bg_profs = osh19_init_bg_profs(params, grid)

bg_profs = struct();

% vertical grid for background moisture profile
zz = transpose(grid.dz:grid.dz:params.params.H-grid.dz);

% Unscaled baroclinic modes
sinz     = transpose(sin(zz * pi / params.H));
re_sinz  = reshape(sinz, [1, 1, length(sinz)]);
sinz_mat = repmat(re_sinz, [params.ny, params.nx, 1]);

sin2z     = transpose(sin(2 * zz * pi / params.H));
re_sin2z  = reshape(sin2z, [1, 1, length(sin2z)]);
sin2z_mat = repmat(re_sin2z, [params.ny, params.nx, 1]);


% Background moisture profile
qbc_bot     = -params.B_vs * params.params.H;
qbc_top     = 0;
atemp       = qbc_bot / (1 - exp(-params.params.H / params.s_tilde));
btemp       = -atemp * exp(-params.H / params.s_tilde);
q_bg_vec_int = atemp * exp(-zz / params.s_tilde) + btemp;
q_bg_vec     = [qbc_bot; q_bg_vec_int; qbc_top];

B_vs_vec     = 1 / (2*grid.dz) * (q_bg_vec(3:end) - q_bg_vec(1:end-2));

% Background moisture damping profiles
tau_vec = params.tau_up + (params.tau_up - params.tau_mid) ...
    * (zz - params.H) / params.H;

D_h_vec = params.D_hUp + (params.D_hUp - params.D_hMid) ...
    * (zz - params.H) / params.H;

D_v_vec = params.D_v + (params.D_v - params.D_v) ...
    * (zz - params.H) / params.H;

bg_profs.sinz     = sinz;
bg_profs.sinz_mat = sinz_mat;

bg_profs.sin2z     = sin2z;
bg_profs.sin2z_mat = sin2z_mat;

bg_profs.q_bg_vec_int = q_bg_vec_int;
bg_profs.q_bg_vec     = q_bg_vec;
bg_profs.B_vs_vec     = B_vs_vec;

bg_profs.tau_vec = tau_vec;
bg_profs.D_h_vec = D_h_vec;
bg_profs.D_v_vec = D_v_vec;

end

