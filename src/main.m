function ierror = main(osh19_params)

ierror = 1;

osh19_params   = osh19_convert_params(osh19_params);
osh19_grid     = osh19_init_grid(osh19_params);
osh19_bg_profs = osh19_init_bg_profs(osh19_params, osh19_grid);
osh19_state    = osh19_init_state(osh19_params, osh19_grid, osh19_bg_profs);

% Output initial condition to file
osh19_output_state(osh19_params, osh19_bg_profs, osh19_state, ...
            0, 0.0);

% Extract parameters to neaten equations
sim_days = osh19_params.sim_days; % "sim_days" is in seconds
out_freq = osh19_params.out_freq;

dt = osh19_grid.dt;

% Set up time-related variables
time         = 0.0;                       % Current time in seconds
out_idx      = 1;                         % Output file number
out_time     = out_idx * out_freq;        % Output file time
days_to_secs = 3600*24;

while time < sim_days 
    osh19_state = osh19_advance_state(osh19_params, osh19_grid, ...
        osh19_bg_profs, osh19_state);
    
    time = time + dt;

    if abs(out_time - time) < dt/2 % Closest time to desired output time
        % check state to ensure validity
        if any(isnan(osh19_state.u))
            out_msg = sprintf(['ERROR: NaNs detected in solution!']);
        end
        
        % write state to file
        osh19_output_state(osh19_params, osh19_bg_profs, osh19_state, ...
            out_idx, time);
        
        % Output message to update user on execution status
        out_msg  = sprintf(['Output file at %.2f days, closest to desired ' ...
            'output time %.2f days.'], ...
            round(time/days_to_secs, 2), ...
            round(out_time/days_to_secs, 2));
        disp(out_msg);
        
        out_idx  = out_idx + 1;
        out_time = out_idx * out_freq;
    end
    
end

end