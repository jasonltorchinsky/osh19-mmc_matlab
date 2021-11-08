function state_out = cmg14_advance_state(params, opers, time, state_in)

state_out = state_in;

% Unpack some common parameters to clean up equations
d_w     = params.d_w;
gamma   = params.gamma;
w_u_hat = params.w_u_hat;

sigma_u = params.sigma_u;
sigma_v = params.sigma_v;
sigma_w = params.sigma_w;

dt      = params.dt;

% IMEX Euler Scheme x_k+1 = (I - dt*A)^(-1) * (x_k + dt*(B + F) + sqrt(dt) W)
% Put into matrices for easier equations

state_vec = [state_out.u_1
    state_out.u_2
    state_out.v
    state_out.w_u];

opers.B = [gamma * state_out.v * state_out.u_1 - state_out.w_u * state_out.u_2
    gamma * state_out.v * state_out.u_2 + state_out.w_u * state_out.u_1
    -gamma * (state_out.u_1^2 + state_out.u_2^2)
    0];

opers.F = [0
    0
    det_force(params, time + dt)
    d_w * w_u_hat];

Sigma = [[sigma_u 0 0 0]
    [0 sigma_u 0 0]
    [0 0 sigma_v 0]
    [0 0 0 sigma_w]];

stoch = normrnd(0, 1, [4,1]);

state_vec = opers.A * (state_vec ...
    + dt * (opers.B + opers.F) ...
    + sqrt(dt) * Sigma * stoch);

% Unpack state vector
state_out.u_1 = state_vec(1);
state_out.u_2 = state_vec(2);
state_out.v   = state_vec(3);
state_out.w_u = state_vec(4);

end

function res = det_force(params, time)

% Unpack some common parameters to clean up equations
d_v = params.d_v;
f_0 = params.f_0;
f_t = params.f_t;
w_f = params.w_f;
phi = params.phi;

res = d_v * (f_0 + f_t * sin(w_f * time + phi)) ...
    + f_t * w_f * cos(w_f * time + phi);

end

