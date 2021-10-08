function state = osh19_prognose_state(params, grid, state_in)

state = state_in;

% Unpack some commonly used variables
nx = params.nx;
ny = params.ny;
nz = params.nz;

g       = params.g;
theta_0 = params.theta_0;

dz = grid.dz;


% Calculate barotropic variables

scale_x = (2 * pi) / (2 * pi * params.P_E);
scale_y = (2 * pi) / (2 * params.P_Y);

k_y_dim_vec  = transpose([0:ny/2 -ny/2+1:-1]) * scale_y;
k_x_dim_vec  = [0:nx/2 -nx/2+1:-1] * scale_x;
[kk_x, kk_y] = meshgrid(k_x_dim_vec, k_y_dim_vec);
wavenum_mat  = -(kk_y.^2+kk_x.^2);
wavenum_mat(1,1) = 10^10;  % THIS NUMBER SHOULD BE IRRELEVANT, AS LONG AS IT'S NONZERO, RIGHT?  
% (THE RHS OF THE PBAR EQUATION SHOULD SATISFY A ZERO MEAN CONDITION SO THAT IT SHOULDN'T MATTER) 

% Calculate *barotropic* winds from barotropic relative vorticity
state.tau_z  = ifft2(fft2(state.zeta_tau)./wavenum_mat);
state.tau_z  = dealias(state.tau_z, 2.0/3.0);
state.u_tau  = -D1(state_in.tau_z,'y',scale_y);
state.v_tau  = D1(state_in.tau_z,'x',scale_x);

% Calculate total horizontal winds at each level
for kk = 2:nz
    state.u(:,:,kk) = state.u_psi(:,:,kk) + state.u_tau;
    state.v(:,:,kk) = state.v_psi(:,:,kk) + state.v_tau;
end 

% Calculate top-level total winds from barotropic winds
state.u(:,:,nz+1) = nz * state.u_tau - squeeze(sum(state.u(:,:,2:nz), 3));
state.v(:,:,nz+1) = nz * state.v_tau - squeeze(sum(state.v(:,:,2:nz), 3));

% Calculate total vertical winds
for kk = 2:nz
    state.w(:,:,kk) = state.w(:,:,kk-1) ...
        - dz * (D1(state.u(:,:,kk),'x',scale_x) ...
        + D1(state.v(:,:,kk),'y',scale_y));
end

% Calculate pressure
for kk = 1:nz-1
    state.p(:,:,2) = state.p(:,:,2) ...
        - g * dz / (nz + 1) * (nz + 1 - kk) ...
        * (1/theta_0 * (state.theta(:,:,kk+1)));
end
for kk = 3:nz+1
    state.p(:,:,kk) = state.p(:,:,kk-1) ...
        + g * dz * (1/theta_0 * (state.theta(:,:,kk-1)));
end

end

