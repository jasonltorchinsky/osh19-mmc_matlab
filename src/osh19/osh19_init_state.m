function state = osh19_init_state(params, grid, bg_profs)

state = struct();

% Unpack some commonly used variables
nx = params.nx;
ny = params.ny;
nz = params.nz;

% Diagnostic variables
state.zeta_tau = zeros([ny, nx]);
state.u_psi    = zeros([ny, nx, nz + 1]);
state.v_psi    = zeros([ny, nx, nz + 1]);
state.theta    = zeros([ny, nx, nz + 1]);
state.q        = zeros([ny, nx, nz + 1]);
% Prognostic variables
state.tau_z    = zeros([ny, nx]);
state.u_tau    = zeros([ny, nx]);
state.v_tau    = zeros([ny, nx]);
state.u        = zeros([ny, nx, nz + 1]);
state.v        = zeros([ny, nx, nz + 1]);
state.w        = zeros([ny, nx, nz + 1]);
state.p        = zeros([ny, nx, nz + 1]);

IC_amps_by_k = [0.7547 0.2760 0.6797 0.6551 0.1626 0.1190 0.4984 0.9597 ...
                0.3404 0.5853 0.2238 0.7513 0.2551 0.5060 0.6991 0.8909 ...
                0.9593 0.5472 0.1386 0.1493];
IC_ps_by_k   = [1.0180 5.2824 1.5977 5.1163 1.5301 5.8387 2.1990 1.2352 ...
                1.5776 3.8707 2.9738 2.2095 5.2203 3.6773 3.4540 5.7629 ...
                1.7960 4.7576 4.7358 2.3904];

for wavenum = params.IC_wavenums
    for mode = params.IC_modes
        lin_state = osh19_calc_lin_state(params, grid, bg_profs, mode, wavenum);
        state = osh19_add_states(1.0, IC_amps_by_k(wavenum), state, lin_state);
    end
end
            
end

