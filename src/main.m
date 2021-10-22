function ierror = main(truth_params, eof_params, dcm_params, mjoo_params)

% Truth is OSH19 system
% EOFs are determined from OSH19 system
% DCM is OSH19 system with modified parameters
% MJOO is modified CMG14 system

ierror = 0; % 0 is no errors

brk_str = '-------------------------------------------------------------------';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Truth simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if truth_params.init_simulation
    disp(brk_str);
    msg = sprintf(['Initializing truth simulation...\n']);
    disp(msg);
    
    % Initialize parameters
    msg = sprintf(['Converting input parameter time units to (s)...\n']);
    disp(msg);
    
    truth_params = osh19_convert_params(truth_params);
    
    osh19_output_params(truth_params);
    
    % Initialize grid
    msg = sprintf(['Initializing grid...\n']);
    disp(msg);
    
    truth_grid = osh19_init_grid(truth_params);
    
    osh19_output_grid(truth_params, truth_grid);
    
    % Initialize background profiles
    msg = sprintf(['Initializing background profiles...\n']);
    disp(msg);
    
    truth_bg_profs = osh19_init_bg_profs(truth_params, truth_grid);
    
    osh19_output_bg_profs(truth_params, truth_bg_profs);
    
    % Initialize state
    msg = sprintf(['Initializing state...\n']);
    disp(msg);
    
    [truth_state, truth_growths_freqs] = osh19_init_state(truth_params, ...
        truth_grid, truth_bg_profs);
    
    osh19_output_growths_freqs(truth_params, truth_growths_freqs);
    
    
    if truth_params.run_simulation
        disp(brk_str);
        msg = sprintf(['Running truth simulation...\n']);
        disp(msg);
        
        % Extract parameters to neaten equations
        sim_days = truth_params.sim_days; % "sim_days" is in seconds
        out_freq = truth_params.out_freq;
        
        dt = truth_grid.dt;
        
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
            1000*max(truth_state.u,[],'all'));
        disp(msg);
        
        osh19_output_state(truth_params, truth_bg_profs, truth_state, ...
            0, 0.0);
        
        msg = sprintf(['Beginning time-stepping...\n']);
        disp(msg);
        
        
        while time < sim_days
            truth_state = osh19_advance_state(truth_params, truth_grid, ...
                truth_bg_profs, truth_state);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % check state to ensure validity
                if any(isnan(truth_state.u))
                    msg = sprintf(['ERROR: NaNs detected in solution!']);
                    disp(msg);
                    
                    ierror = 1;
                    
                    return
                end
                
                % write state to file
                osh19_output_state(truth_params, truth_bg_profs, ...
                    truth_state, out_idx, time);
                
                % Output message to update user on execution status
                msg  = sprintf(['Day %.2f of %.2f.\n'...
                    '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    1000*max(truth_state.u,[],'all'));
                disp(msg);
                
                out_idx  = out_idx + 1;
                out_time = out_idx * out_freq;
                
            end
            
        end
        
    end
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% EOF calculation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if eof_params.calc_eofs
    disp(brk_str);
    msg = sprintf(['Calculating EOFs...\n']);
    disp(msg);
    
    eofs = eof_calc_eof(eof_params);
    
    eof_output_eofs(eof_params, eofs);
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DCM simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if dcm_params.init_simulation
    disp(brk_str);
    msg = sprintf(['Initializing DCM simulation...\n']);
    disp(msg);
    
    % Initialize parameters
    msg = sprintf(['Converting input parameter time units to (s)...\n']);
    disp(msg);
    
    dcm_params = osh19_convert_params(dcm_params);
    
    osh19_output_params(dcm_params);
    
    % Initialize grid
    msg = sprintf(['Initializing grid...\n']);
    disp(msg);
    
    dcm_grid = osh19_init_grid(dcm_params);
    
    osh19_output_grid(dcm_params, dcm_grid);
    
    % Initialize background profiles
    msg = sprintf(['Initializing background profiles...\n']);
    disp(msg);
    
    dcm_bg_profs = osh19_init_bg_profs(dcm_params, dcm_grid);
    
    osh19_output_bg_profs(dcm_params, dcm_bg_profs);
    
    % Initialize state
    msg = sprintf(['Initializing state...\n']);
    disp(msg);
    
    [dcm_state, dcm_growths_freqs] = osh19_init_state(dcm_params, ...
        dcm_grid, dcm_bg_profs);
    
    osh19_output_growths_freqs(dcm_params, dcm_growths_freqs);
    
    
    if dcm_params.run_simulation
        disp(brk_str);
        msg = sprintf(['Running DCM simulation...\n']);
        disp(msg);
        
        % Extract parameters to neaten equations
        sim_days = dcm_params.sim_days; % "sim_days" is in seconds
        out_freq = dcm_params.out_freq;
        
        dt = dcm_grid.dt;
        
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
            1000*max(dcm_state.u,[],'all'));
        disp(msg);
        
        osh19_output_state(dcm_params, dcm_bg_profs, dcm_state, ...
            0, 0.0);
        
        msg = sprintf(['Beginning time-stepping...\n']);
        disp(msg);
        
        while time < sim_days
            dcm_state = osh19_advance_state(dcm_params, dcm_grid, ...
                dcm_bg_profs, dcm_state);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % check state to ensure validity
                if any(isnan(dcm_state.u))
                    msg = sprintf(['ERROR: NaNs detected in solution!']);
                    disp(msg);
                    
                    ierror = 1;
                    
                    return
                end
                
                % write state to file
                osh19_output_state(dcm_params, dcm_bg_profs, ...
                    dcm_state, out_idx, time);
                
                % Output message to update user on execution status
                msg  = sprintf(['Day %.2f of %.2f.\n'...
                    '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    1000*max(dcm_state.u,[],'all'));
                disp(msg);
                
                out_idx  = out_idx + 1;
                out_time = out_idx * out_freq;
                
            end
            
        end
        
    end
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% MJOO simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if mjoo_params.init_simulation
    disp(brk_str);
    msg = sprintf(['Initializing MJO-only simulation...\n']);
    disp(msg);
    
    % Initialize parameters
    msg = sprintf(['Converting input parameter time units to (s)...\n']);
    disp(msg);
    
    mjoo_params = cmg14_convert_params(mjoo_params);
    
    cmg14_output_params(mjoo_params);

    % Initialize time-stepping operators
    msg = sprintf(['Initializing background profiles...\n']);
    disp(msg);
    
    mjoo_opers = cmg14_init_opers(mjoo_params);
    
    % Initialize state
    msg = sprintf(['Initializing state...\n']);
    disp(msg);
    
    mjoo_state = cmg14_init_state(mjoo_params);
    
    if mjoo_params.run_simulation
        disp(brk_str);
        msg = sprintf(['Running MJO-only simulation...\n']);
        disp(msg);
        
        % Extract parameters to neaten equations
        sim_days = mjoo_params.sim_days; % "sim_days" is in seconds
        out_freq = mjoo_params.out_freq;
        
        dt = mjoo_params.dt;
        
        % Set up time-related variables
        time         = 0.0;                 % Current time in seconds
        out_idx      = 1;                   % Output file number
        out_time     = out_idx * out_freq;  % Output file time
        days_to_secs = 3600*24;
        
        % Output initial state
        
        msg  = sprintf(['Day %.2f of %.2f.\n'...
            'State: u_1 %.4f' ...
            '       u_2 %.4f' ...
            '       v   %.4f' ...
            '       w_u %.4f'], ...
            round(0.0/days_to_secs, 2), ...
            round(0/days_to_secs, 2), ...
            round(mjoo_state.u_1, 4), ...
            round(mjoo_state.u_2, 4), ...
            round(mjoo_state.v, 4), ...
            round(mjoo_state.w_u, 4));
        disp(msg);
        
        cmg14_output_state(mjoo_params, mjoo_state, 0, 0.0);
        
        msg = sprintf(['Beginning time-stepping...\n']);
        disp(msg);
        
        while time < sim_days
            mjoo_state = cmg14_advance_state(mjoo_params, mjoo_opers, time, ...
                mjoo_state);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % Check state to ensure validity
                
                % write state to file
                cmg14_output_state(mjoo_params, mjoo_state, out_idx, time);
                
                % Output message to update user on execution status
                msg  = sprintf(['Day %.2f of %.2f.\n'...
                    'State: u_1 %.4f' ...
                    '       u_2 %.4f' ...
                    '       v   %.4f' ...
                    '       w_u %.4f'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    round(mjoo_state.u_1, 4), ...
                    round(mjoo_state.u_2, 4), ...
                    round(mjoo_state.v, 4), ...
                    round(mjoo_state.w_u, 4));
                disp(msg);
                
                out_idx  = out_idx + 1;
                out_time = out_idx * out_freq;
                
            end
            
        end
        
    end
    
end


end