function ierror = main(osh19_params)

ierror = 1;

% Initialize parameters
msg = sprintf(['Converting input parameters to time units (s)...\n']);
disp(msg);

osh19_params   = osh19_convert_params(osh19_params);

% Initialize grid
msg = sprintf(['Initializing grid...\n']);
disp(msg);

osh19_grid     = osh19_init_grid(osh19_params);

msg = sprintf(['Writing grid to file...\n']);
disp(msg);

osh19_output_grid(osh19_params, osh19_grid);

% Initialize bacgkround profiles
msg = sprintf(['Initializing background profiles...\n']);
disp(msg);

osh19_bg_profs = osh19_init_bg_profs(osh19_params, osh19_grid);

pause(inf);

% Initialize state
msg = sprintf(['Initializing state...\n']);
disp(msg);

osh19_state    = osh19_init_state(osh19_params, osh19_grid, osh19_bg_profs);

msg = sprintf(['Writing initial condition to file...\n']);
disp(msg);

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
            msg = sprintf(['ERROR: NaNs detected in solution!']);
            disp(msg);
        end
        
        % write state to file
        osh19_output_state(osh19_params, osh19_bg_profs, osh19_state, ...
            out_idx, time);
        
        % Output message to update user on execution status
        msg  = sprintf(['Output file at %.2f days,\n'...
            '   closest to desired output time %.2f days.'], ...
            round(time/days_to_secs, 2), ...
            round(out_time/days_to_secs, 2));
        disp(msg);
        
        out_idx  = out_idx + 1;
        out_time = out_idx * out_freq;
    end
    
end

end