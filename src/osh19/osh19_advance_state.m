function state_out = osh19_advance_state(params, grid, bg_profs, state_in)

% Implements RK4
% Y1 = U^n
% Y2 = U^n + dt/2 * fY1
% Y3 = U^n + dt/2 * fY2
% Y4 = U^n + dt * fY3
% U^(n+1) = U^n + dt/6 * fY1 + dt/3 * fY2 + dt/3 * fY3 + dt/6 * fY4

dt = grid.dt;

% Dealias the input state
stage_temp = osh19_dealias_state(state_in, 2.0/3.0); % Holds intermediate states

% Perform the RK4 update
fY1        = osh19_calc_RK4_stage(params, grid, bg_profs, stage_temp);
stage_temp = osh19_add_states(1.0, dt/2.0, state_in, fY1); %Y2

stage_temp = osh19_dealias_state(stage_temp, 2.0/3.0);
fY2        = osh19_calc_RK4_stage(params, grid, bg_profs, stage_temp);
stage_temp = osh19_add_states(1.0, dt/2.0, state_in, fY2); %Y3

stage_temp = osh19_dealias_state(stage_temp, 2.0/3.0);
fY3        = osh19_calc_RK4_stage(params, grid, bg_profs, stage_temp);
stage_temp = osh19_add_states(1.0, dt, state_in, fY3); %Y4

stage_temp = osh19_dealias_state(stage_temp, 2.0/3.0);
fY4        = osh19_calc_RK4_stage(params, grid, bg_profs, stage_temp);

state_out = osh19_add_states(1.0, dt/6.0, state_in,  fY1);
state_out = osh19_add_states(1.0, dt/3.0, state_out, fY2);
state_out = osh19_add_states(1.0, dt/3.0, state_out, fY3);
state_out = osh19_add_states(1.0, dt/6.0, state_out, fY4);


end

