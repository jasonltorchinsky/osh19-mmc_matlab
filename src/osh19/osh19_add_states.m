function state = osh19_add_states(a, b, state_a, state_b)

state = struct();
state.zeta_tau = a * state_a.zeta_tau + b * state_b.zeta_tau;
state.u_clin   = a * state_a.u_clin   + b * state_b.u_clin;
state.v_clin   = a * state_a.v_clin   + b * state_b.v_clin;
state.theta    = a * state_a.theta    + b * state_b.theta;
state.q        = a * state_a.q        + b * state_b.q;
state.tau_z    = a * state_a.tau_z    + b * state_b.tau_z;
state.u_tau    = a * state_a.u_tau    + b * state_b.u_tau;
state.v_tau    = a * state_a.v_tau    + b * state_b.v_tau;
state.u        = a * state_a.u        + b * state_b.u;
state.v        = a * state_a.v        + b * state_b.v;
state.w        = a * state_a.w        + b * state_b.w;
state.p        = a * state_a.p        + b * state_b.p;

end

