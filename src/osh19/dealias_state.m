function state_out = dealias_state(state_in, fraction)

state_out = state_in;

state_out.zeta_tau = dealias(state_out.zeta_tau, fraction);
state_out.u_psi    = dealias(state_out.u_psi,    fraction);
state_out.v_psi    = dealias(state_out.v_psi,    fraction);
state_out.theta    = dealias(state_out.theta,    fraction);
state_out.q        = dealias(state_out.q,        fraction);

state_out.tau_z = dealias(state_out.tau_z, fraction);
state_out.u_tau = dealias(state_out.u_tau, fraction);
state_out.v_tau = dealias(state_out.v_tau, fraction);
state_out.u     = dealias(state_out.u,     fraction);
state_out.v     = dealias(state_out.v,     fraction);
state_out.w     = dealias(state_out.w,     fraction);
state_out.p     = dealias(state_out.p,     fraction);

end

