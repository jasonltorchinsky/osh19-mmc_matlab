function ierror = main(truth_params, eof_params, dcm_params, mjoo_params, ...
    ens_params)

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
    fprintf('Initializing truth simulation...\n');
    
    % Initialize parameters
    fprintf('Converting input parameter time units to (s)...\n');
    
    truth_params = osh19_convert_params(truth_params);
    
    osh19_output_params(truth_params);
    
    % Initialize grid
    fprintf('Initializing grid...\n');
    
    truth_grid = osh19_init_grid(truth_params);
    
    osh19_output_grid(truth_params, truth_grid);
    
    % Initialize background profiles
    fprintf('Initializing background profiles...\n');
    
    
    truth_bg_profs = osh19_init_bg_profs(truth_params, truth_grid);
    
    osh19_output_bg_profs(truth_params, truth_bg_profs);
    
    % Initialize state
    fprintf('Initializing state...\n');
    
    
    [truth_state, truth_growths_freqs] = osh19_init_state(truth_params, ...
        truth_grid, truth_bg_profs);
    
    osh19_output_growths_freqs(truth_params, truth_growths_freqs);
    
    
    if truth_params.run_simulation
        disp(brk_str);
        fprintf('Running truth simulation...\n');
        
        
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
        
        fprintf(['Day %.2f of %.2f.\n'...
            '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
            round(0.0/days_to_secs, 2), ...
            round(sim_days/days_to_secs, 2), ...
            1000*max(truth_state.u,[],'all'));
        
        
        osh19_output_state(truth_params, truth_bg_profs, truth_state, ...
            0, 0.0);
        
        fprintf('Beginning time-stepping...\n');
        
        
        
        while time < sim_days
            truth_state = osh19_advance_state(truth_params, truth_grid, ...
                truth_bg_profs, truth_state);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % check state to ensure validity
                if any(isnan(truth_state.u))
                    fprintf('ERROR: NaNs detected in solution!');
                    
                    
                    ierror = 1;
                    
                    return
                end
                
                % write state to file
                osh19_output_state(truth_params, truth_bg_profs, ...
                    truth_state, out_idx, time);
                
                % Output message to update user on execution status
                fprintf(['Day %.2f of %.2f.\n'...
                    '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    1000*max(truth_state.u,[],'all'));
                
                
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
    fprintf('Calculating EOFs...\n');
    
    
    eofs = eof_calc_eof(eof_params);
    
    eof_output_eofs(eof_params, eofs);
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DCM simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if dcm_params.init_simulation
    disp(brk_str);
    fprintf('Initializing DCM simulation...\n');
    
    
    % Initialize parameters
    fprintf('Converting input parameter time units to (s)...\n');
    
    
    dcm_params = osh19_convert_params(dcm_params);
    
    osh19_output_params(dcm_params);
    
    % Initialize grid
    fprintf('Initializing grid...\n');
    
    
    dcm_grid = osh19_init_grid(dcm_params);
    
    osh19_output_grid(dcm_params, dcm_grid);
    
    % Initialize background profiles
    fprintf('Initializing background profiles...\n');
    
    
    dcm_bg_profs = osh19_init_bg_profs(dcm_params, dcm_grid);
    
    osh19_output_bg_profs(dcm_params, dcm_bg_profs);
    
    % Initialize state
    fprintf('Initializing state...\n');
    
    
    [dcm_state, dcm_growths_freqs] = osh19_init_state(dcm_params, ...
        dcm_grid, dcm_bg_profs);
    
    osh19_output_growths_freqs(dcm_params, dcm_growths_freqs);
    
    
    if dcm_params.run_simulation
        disp(brk_str);
        fprintf('Running DCM simulation...\n');
        
        
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
        
        fprintf(['Day %.2f of %.2f.\n'...
            '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
            round(0.0/days_to_secs, 2), ...
            round(sim_days/days_to_secs, 2), ...
            1000*max(dcm_state.u,[],'all'));
        
        
        osh19_output_state(dcm_params, dcm_bg_profs, dcm_state, ...
            0, 0.0);
        
        fprintf('Beginning time-stepping...\n');
        
        
        while time < sim_days
            dcm_state = osh19_advance_state(dcm_params, dcm_grid, ...
                dcm_bg_profs, dcm_state);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % check state to ensure validity
                if any(isnan(dcm_state.u))
                    fprintf('ERROR: NaNs detected in solution!');
                    
                    
                    ierror = 1;
                    
                    return
                end
                
                % write state to file
                osh19_output_state(dcm_params, dcm_bg_profs, ...
                    dcm_state, out_idx, time);
                
                % Output message to update user on execution status
                fprintf(['Day %.2f of %.2f.\n'...
                    '   Max zonal wind speed: %.4f m s^(-1).\n'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    1000*max(dcm_state.u,[],'all'));
                
                
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
    fprintf('Initializing MJO-only simulation...\n');
    
    % Initialize parameters
    fprintf('Converting input parameter time units to (s)...\n');
    
    mjoo_params = cmg14_convert_params(mjoo_params);
    
    if (mjoo_params.dt > mjoo_params.out_freq)
        fprintf(['WARNING: Desired output frequency is less'...
            ' than step size. Adjust output frequency...\n']);
        
        mjoo_params.out_freq = mjoo_params.dt;
    end
    
    cmg14_output_params(mjoo_params);

    % Initialize time-stepping operators
    fprintf('Initializing time-stepping operators...\n');
    
    mjoo_opers = cmg14_init_opers(mjoo_params);
    
    % Initialize state
    fprintf('Initializing state...\n');
    
    
    mjoo_state = cmg14_init_state(mjoo_params);
    
    if mjoo_params.run_simulation
        disp(brk_str);
        fprintf('Running MJO-only simulation...\n');
        
        
        % Extract parameters to neaten equations
        sim_days = mjoo_params.sim_days; % "sim_days" is in seconds
        out_freq = mjoo_params.out_freq; % out_freq is in seconds
        
        dt = mjoo_params.dt;
        

        % Set up time-related variables
        time         = 0.0;                 % Current time in seconds
        out_idx      = 1;                   % Output file number
        out_time     = out_idx * out_freq;  % Output file time
        days_to_secs = 3600*24;
        
        % Output initial state

        fprintf(['Day %.2f of %.2f.\n'...
            'State: u_1 %.4f' ...
            '       u_2 %.4f' ...
            '       v   %.4e' ...
            '       w_u %.4e\n'], ...
            round(0.0/days_to_secs, 2), ...
            round(0/days_to_secs, 2), ...
            round(mjoo_state.u_1, 4), ...
            round(mjoo_state.u_2, 4), ...
            round(mjoo_state.v, 4, 'significant'), ...
            round(mjoo_state.w_u, 4, 'significant'));
        
        
        cmg14_output_state(mjoo_params, mjoo_state, 0, 0.0);
        
        fprintf('Beginning time-stepping...\n');
        
        
        while time < sim_days
            mjoo_state = cmg14_advance_state(mjoo_params, mjoo_opers, time, ...
                mjoo_state);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % Check state to ensure validity
                
                % write state to file
                cmg14_output_state(mjoo_params, mjoo_state, out_idx, time);
                
                % Output message to update user on execution status
                fprintf(['Day %.2f of %.2f.\n'...
                    'State: u_1 %.4f' ...
                    '       u_2 %.4f' ...
                    '       v   %.4e' ...
                    '       w_u %.4e\n'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    round(mjoo_state.u_1, 4), ...
                    round(mjoo_state.u_2, 4), ...
                    round(mjoo_state.v, 4, 'significant'), ...
                    round(mjoo_state.w_u, 4, 'significant'));
                
                
                out_idx  = out_idx + 1;
                out_time = out_idx * out_freq;
                
            end
            
        end
        
    end
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Multi-Model Ensemble
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ens_params.init_simulation
    disp(brk_str);
    fprintf('Initializing ensemble simulation...\n');
    
    
    % Unpack DCM, MJOO parameters. From here, DO NOT TOUCH DCM, MJOO PARAMS
    % WITHIN ENS_PARAMS
    ens_dcm_params = ens_params.dcm_params;
    ens_mjoo_params = ens_params.mjoo_params;
    
    fprintf('Converting input parameter time units to (s)...\n');

    ens_dcm_params = osh19_convert_params(ens_dcm_params);
    osh19_output_params(ens_dcm_params);
    
    ens_mjoo_params = cmg14_convert_params(ens_mjoo_params);
    cmg14_output_params(ens_mjoo_params);
    
    % Initialize DCM component
    fprintf('Initializing DCM component...\n');
    
    ens_dcm_grid = osh19_init_grid(ens_dcm_params);
    osh19_output_grid(ens_dcm_params, ens_dcm_grid);
    
    ens_dcm_bg_profs = osh19_init_bg_profs(ens_dcm_params, ens_dcm_grid);
    osh19_output_bg_profs(ens_dcm_params, ens_dcm_bg_profs);
    
    [ens_dcm_state, ens_dcm_growths_freqs] = osh19_init_state(ens_dcm_params, ...
        ens_dcm_grid, ens_dcm_bg_profs);
    osh19_output_growths_freqs(ens_dcm_params, ens_dcm_growths_freqs);
    
    % Initialize MJOO component
    fprintf('Initializing MJOO component...\n');
    
    % Copy dt, sim_days, out_freq from DCM component
    ens_mjoo_params.dt = ens_dcm_grid.dt;
    ens_mjoo_params.sim_days = ens_dcm_params.sim_days;
    ens_mjoo_params.out_freq = ens_dcm_params.out_freq;
    
    ens_mjoo_opers = cmg14_init_opers(ens_mjoo_params);
    ens_mjoo_state = cmg14_init_state(ens_mjoo_params);
    
    % Read-in EOFs for communication information
    ens_eofs = ens_read_eofs(ens_params);

    
    if ens_params.run_simulation
        disp(brk_str);
        fprintf('Running ensemble simulation...\n');
        
        
        % Extract parameters to neaten equations
        % These are the same between both components
        dt = ens_dcm_grid.dt;
        sim_days = ens_dcm_params.sim_days; % sim_days is in seconds
        out_freq = ens_dcm_params.out_freq; % out_freq is in seconds
        
        
        % Set up time-related variables
        time         = 0.0;                 % Current time in seconds
        out_idx      = 1;                   % Output file number
        out_time     = out_idx * out_freq;  % Output file time
        days_to_secs = 3600*24;
        
        % Output initial state
        fprintf(['Day %.2f of %.2f.\n'...
            '   DCM:  Max zonal wind speed: %.4f m s^(-1).\n'...
            '   MJOO: u_1 %.4f' ...
            '         u_2 %.4f' ...
            '         v   %.4e' ...
            '         w_u %.4e\n'], ...
            round(0.0/days_to_secs, 2), ...
            round(sim_days/days_to_secs, 2), ...
            1000*max(ens_dcm_state.u,[],'all'), ...
            round(ens_mjoo_state.u_1, 4), ...
            round(ens_mjoo_state.u_2, 4), ...
            round(ens_mjoo_state.v, 4, 'significant'), ...
            round(ens_mjoo_state.w_u, 4, 'significant'));
        
        
        osh19_output_state(ens_dcm_params, ens_dcm_bg_profs, ens_dcm_state, ...
            0, 0.0);
        cmg14_output_state(ens_mjoo_params, ens_mjoo_state, 0, 0.0);
        
        fprintf('Beginning time-stepping...\n');
        
        
        while time < sim_days
            ens_dcm_state = osh19_advance_state(ens_dcm_params, ens_dcm_grid, ...
                ens_dcm_bg_profs, ens_dcm_state);
            ens_mjoo_state = cmg14_advance_state(ens_mjoo_params, ens_mjoo_opers, time, ...
                ens_mjoo_state);
            
            [ens_dcm_state, ens_mjoo_state] = ens_comm(dcm_params, ...
                ens_dcm_grid, ens_dcm_state, ens_mjoo_state, ...
                ens_params, ens_eofs);
            
            time = time + dt;
            
            if abs(out_time - time) < dt/2 % Closest time to desired output time
                % check state to ensure validity
                if any(isnan(ens_dcm_state.u))
                    fprintf('ERROR: NaNs detected in solution!');
                    
                    ierror = 1;
                    
                    return
                end
                
                % Write state to file
                osh19_output_state(ens_dcm_params, ens_dcm_bg_profs, ...
                    ens_dcm_state, out_idx, time);
                cmg14_output_state(ens_mjoo_params, ens_mjoo_state, out_idx, time);
                
                % Output message to update user on execution status
                fprintf(['Day %.2f of %.2f.\n'...
                    '   DCM:  Max zonal wind speed: %.4f m s^(-1).\n'...
                    '   MJOO: u_1 %.4f' ...
                    '         u_2 %.4f' ...
                    '         v   %.4e' ...
                    '         w_u %.4e\n'], ...
                    round(time/days_to_secs, 2), ...
                    round(sim_days/days_to_secs, 2), ...
                    1000*max(ens_dcm_state.u,[],'all'), ...
                    round(ens_mjoo_state.u_1, 4), ...
                    round(ens_mjoo_state.u_2, 4), ...
                    round(ens_mjoo_state.v, 4, 'significant'), ...
                    round(ens_mjoo_state.w_u, 4, 'significant'));
                
                out_idx  = out_idx + 1;
                out_time = out_idx * out_freq;
                
            end
            
        end
        
    end
    
end


end