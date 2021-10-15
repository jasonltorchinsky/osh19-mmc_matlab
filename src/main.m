function ierror = main(truth_params, eof_params)

% Truth is OSH19 system
% EOFs are determined from OSH19 system
% DCM is OSH19 system with modified parameters
% MJOO is modified CMG14 system

ierror = 0; % 0 is no errors

% Truth simulation
if truth_params.init_simulation
    
    % Initialize parameters
    msg = sprintf(['Converting input parameter time units to (s)...\n']);
    disp(msg);
    
    truth_params = osh19_convert_params(truth_params);
    
    osh19_output_params(truth_params);
    
    % Initialize grid
    msg = sprintf(['Initializing grid...\n']);
    disp(msg);
    
    osh19_grid = osh19_init_grid(truth_params);
    
    osh19_output_grid(truth_params, osh19_grid);
    
    % Initialize bacgkround profiles
    msg = sprintf(['Initializing background profiles...\n']);
    disp(msg);
    
    osh19_bg_profs = osh19_init_bg_profs(truth_params, osh19_grid);
    
    osh19_output_bg_profs(truth_params, osh19_bg_profs);
    
    % Initialize state
    msg = sprintf(['Initializing state...\n']);
    disp(msg);
    
    [osh19_state, osh19_growths_freqs] = osh19_init_state(truth_params, ...
        osh19_grid, osh19_bg_profs);
    
    osh19_output_growths_freqs(truth_params, osh19_growths_freqs);
    
    
    if truth_params.run_simulation
        
        % Extract parameters to neaten equations
        sim_days = truth_params.sim_days; % "sim_days" is in seconds
        out_freq = truth_params.out_freq;
        
        dt = osh19_grid.dt;
        
        % Set up time-related variables
        time         = 0.0;                 % Current time in seconds
        out_idx      = 1;                   % Output file number
        out_time     = out_idx * out_freq;  % Output file time
        days_to_secs = 3600*24;
        
        % Output initial state
        
        msg  = sprintf(['Day %.2f of %.2f.\n'...
            '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
            round(0.0/days_to_secs, 2), ...
            round(sim_days/days_to_secs, 2), ...
            1000*max(osh19_state.u,[],'all'));
        disp(msg);
        
        osh19_output_state(truth_params, osh19_bg_profs, osh19_state, ...
            0, 0.0);
        
        msg = sprintf(['Beginning time-stepping...\n']);
        disp(msg);
        
        while time < sim_days
            osh19_state = osh19_advance_state(truth_params, osh19_grid, ...
                osh19_bg_profs, osh19_state);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % check state to ensure validity
                if any(isnan(osh19_state.u))
                    msg = sprintf(['ERROR: NaNs detected in solution!']);
                    disp(msg);
                    
                    ierror = 1;
                    
                    return
                end
                
                % write state to file
                osh19_output_state(truth_params, osh19_bg_profs, ...
                    osh19_state, out_idx, time);
                
                % Output message to update user on execution status
                msg  = sprintf(['Day %.2f of %.2f.\n'...
                    '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    1000*max(osh19_state.u,[],'all'));
                disp(msg);
                
                out_idx  = out_idx + 1;
                out_time = out_idx * out_freq;
                
            end
            
        end
        
    end
    
end

% EOF calculation
if eof_params.calc_eofs
    
    msg = sprintf(['Calculating EOFs...\n']);
    disp(msg);
    
    eofs = eof_calc_eof(eof_params);
    
    eof_output_eofs(eof_params, eofs);
    
end

% Run ensemble

end