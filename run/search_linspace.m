% Search a range of parameter values for the leading growthrates

out_path = 'output';
exp_name = 'param_search';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Base parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

params = struct();

params.nx = 32; % Number of points zonally
params.ny = 32; % Number of points meridionally
params.nz = 8;  % Number of non-trivial levels for u, v, p
params.H = 16;  % Height of troposphere (km)
params.P_Y = 6000; % Distance from equator to channel wall (km)
params.P_E = 40000/(2*pi); % Radius of planet (km)
params.beta = 2.3 * 10^(-8); % Meridional variation of Coriolis parameter (km^(-1) s^(-1))
params.g = 9.8 * 10^(-3); % Acceleration due to gravity (km s^(-2))
params.c_p = 1000; % Specific heat of dry air at constant pressure (J kg^(-1) K^(-1))
params.L_v = 2.5 * 10^(6); % Latent heat of vaporization (J kg^(-1))
params.B = 3; % Background potential temperature vertical gradient (K km^(-1))
params.theta_0 = 300; % Reference potential temperature (K)
params.tau_u = 16; % Wind damping time-scale (d)
params.tau_theta = 16; % Potential temperature damping time-scale (d)
params.tau_up = 1; % Moisture damping time-scale in upper-troposphere (d)
params.tau_mid = 4/24; % Moisture damping time-scale in mid-troposphere (d)
params.B_vs = -1.2 * 10^(-3); % Mean vertical q_bg gradient (kg kg^(-1) km^(-1))
params.a = 0.2; % 1--pole-to-equator q_bg ratio
params.L_tilde = 3000; % q_bg meridional decay length scale (km)
params.s_tilde = 12; % q_bg vertical decay length scale (km)
params.D_hUp = 60.8; % Horizontal q diffusion in upper-troposphere (km^(2) s^(-1))
params.D_hMid = 7.6; % Horizontal q diffusion in mid-troposphere (km^(2) s^(-1))
params.D_v = 2 * 10^(-4); % Vertical q diffusion (km^(2) s^(-1))
params.IC_type = 1; % Initial condition type:
                          % 1 = get modes, wavenumbers from linear solution
                          % 2 = load state from file
params.IC_modes = [1]; % Modes to use for initial condition (IC)
params.IC_wavenums = 1:1:7; % Zonal wavenumber to use for IC
params.IC_amp = 1; % Amplification factor for IC
params.clin_conv_adj = 1; % Options for baroclinic modes in IC
params.sim_days = sim_days; % Number of days to simulate (d)
params.out_freq = 1; % How often to output data (d)

params.out_path = out_path;
params.exp_name = exp_name;
params.component_name = 'truth';

params.init_simulation = true;
params.run_simulation = false;
params.create_plots = false;

params_base = params;

n_tests = 9;
B_vss         = linspace(-1.5 * 10^(-3), -1.0 * 10^(-3), n_tests);
D_h_Multiples = linspace(0.05, 2.0, n_tests);
D_v_Multiples = linspace(0.05, 2.0, n_tests);

n_modes = size(params.IC_modes, 2);
n_wavenums = size(params.IC_wavenums, 2);
growths = zeros([n_wavenums, n_tests, n_tests, n_tests]);
freqs   = zeros([n_wavenums, n_tests, n_tests, n_tests]);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Solve the linearized system a whole bunch of times
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for B_vs_idx = 1:n_tests
    B_vs = B_vss(B_vs_idx);
    
    for D_h_Multiple_idx = 1:n_tests        
        D_h_Multiple = D_h_Multiples(D_h_Multiple_idx);
        
        for D_v_Multiple_idx = 1:n_tests
            D_v_Multiple = D_v_Multiples(D_v_Multiple_idx);
        
            test_num = (B_vs_idx-1)*n_tests^2 ...
                + (D_h_Multiple_idx-1)*n_tests ...
                + D_v_Multiple_idx;
            fprintf('Test %d of %d.\n', test_num, n_tests^3);
        
            params = params_base;
            params.B_vs   = B_vs;
            params.D_hUp  = D_h_Multiple * params.D_hUp;
            params.D_hMid = D_h_Multiple * params.D_hMid;
            params.D_v    = D_v_Multiple * params.D_v;
            
            params = osh19_convert_params(params);
            grid = osh19_init_grid(params);
            bg_profs = osh19_init_bg_profs(params, grid);
            
            tic;
            for mode_idx = 1:n_modes
                mode = params.IC_modes(mode_idx);
                
                for wavenum_idx = 1:n_wavenums
                    wavenum = params.IC_wavenums(wavenum_idx);
                    
                    [~, growth_rates, freq_rates] = osh19_calc_lin_state(params, ...
                        grid, bg_profs, mode, wavenum);
                growths(wavenum_idx, B_vs_idx, D_h_Multiple_idx, D_v_Multiple_idx) ...
                    = growth_rates(1);
                freqs(wavenum_idx, B_vs_idx, D_h_Multiple_idx, D_v_Multiple_idx) ...
                    = freq_rates(1);
                end
            end
            elapsed_time = toc;
        
            fprintf('Elapsed time (s): %.4f\n', elapsed_time);
            if test_num == 1
                fprintf('Estimated execution time (min): %.4f\n', ...
                    elapsed_time * n_tests^3 / 60);
            end
            
            fprintf('\n');
        
        end
    end
end

levels = [-0.25 -0.2 -0.15 -0.1 -0.05 0 0.25 0.5 0.75 1 1.5 2 3 4];

for wavenum_idx = 1:n_wavenums
    for D_v_Multiple_idx = 1:n_tests
        D_v_Multiple = D_v_Multiples(D_v_Multiple_idx);

        clf;
        wavenum = params.IC_wavenums(wavenum_idx);
    
        contourf(D_h_Multiples, B_vss, ...
            squeeze(growths(wavenum, :, :, D_v_Multiple_idx)), ...
            levels, ShowText = 'on');
        colorbar;
        colormap cool;
        
        title(sprintf('Growthrates for Zonal Wavenumber %d', wavenum));
        txt = texlabel('D_v');
        subtitle_str = strcat([txt sprintf(' Multiple %.2f', D_v_Multiple)]);
        subtitle(subtitle_str);
    
        xlabel(texlabel('D_{h} Multiple'));
        ylabel(texlabel('B_{vs} (kg kg^{-1} km^{-1})'));
    
        hold on;
        plot(D_h_Multiples, params_base.B_vs*ones(size(D_h_Multiples)), 'k-');
        plot(ones(size(B_vss)), B_vss, 'k-');
        hold off;
        
        out_path       = params.out_path;
        exp_path       = fullfile(out_path, params.exp_name);
    
        if ~(isfolder(exp_path))
            mkdir(exp_path);
        end
    
        addpath(out_path);
        addpath(exp_path);
    
        set(gcf, 'PaperUnits', 'inches');
        figWidth  = 7.5; % Figure width in inches.
        figHeight = 6; % Figure height in inches.
        set(gcf,...
            'PaperPosition', [0, 0, figWidth, figHeight],...
            'PaperSize', [figWidth, figHeight])
    
        file_name = sprintf('growthrates_Bvs_Dh_%d_%d.pdf', ...
            wavenum, D_v_Multiple_idx);
        plot_file = fullfile(exp_path, file_name);
        print(plot_file, '-dpdf', '-painters');
    end
end

for wavenum_idx = 1:n_wavenums
    for D_h_Multiple_idx = 1:n_tests
        D_h_Multiple = D_h_Multiples(D_h_Multiple_idx);

        clf;
        wavenum = params.IC_wavenums(wavenum_idx);
    
        contourf(D_v_Multiples, B_vss, ...
            squeeze(growths(wavenum, :, D_h_Multiple_idx, :)), ...
            levels, ShowText = 'on');
        colorbar;
        colormap cool;
        
        title(sprintf('Growthrates for Zonal Wavenumber %d', wavenum));
        txt = texlabel('D_h');
        subtitle_str = strcat([txt sprintf(' Multiple %.2f', D_h_Multiple)]);
        subtitle(subtitle_str);
    
        xlabel(texlabel('D_{v} Multiple'));
        ylabel(texlabel('B_{vs} (kg kg^{-1} km^{-1})'));
    
        hold on;
        plot(D_v_Multiples, params_base.B_vs*ones(size(D_v_Multiples)), 'k-');
        plot(ones(size(B_vss)), B_vss, 'k-');
        hold off;
        
        out_path       = params.out_path;
        exp_path       = fullfile(out_path, params.exp_name);
    
        if ~(isfolder(exp_path))
            mkdir(exp_path);
        end
    
        addpath(out_path);
        addpath(exp_path);
    
        set(gcf, 'PaperUnits', 'inches');
        figWidth  = 7.5; % Figure width in inches.
        figHeight = 6; % Figure height in inches.
        set(gcf,...
            'PaperPosition', [0, 0, figWidth, figHeight],...
            'PaperSize', [figWidth, figHeight])
    
        file_name = sprintf('growthrates_Bvs_Dv_%d_%d.pdf', ...
            wavenum, D_h_Multiple_idx);
        plot_file = fullfile(exp_path, file_name);
        print(plot_file, '-dpdf', '-painters');
    end
end

for wavenum_idx = 1:n_wavenums
    for B_vs_idx = 1:n_tests
        B_vs = B_vss(B_vs_idx);

        clf;
        wavenum = params.IC_wavenums(wavenum_idx);
    
        contourf(D_h_Multiples, D_v_Multiples, ...
            squeeze(growths(wavenum, B_vs_idx, :, :)), ...
            levels, ShowText = 'on');
        colorbar;
        colormap cool;
        
        title(sprintf('Growthrates for Zonal Wavenumber %d', wavenum));
        txt = texlabel('B_{vs}');
        subtitle_str = strcat([txt sprintf(' %.2e', B_vs)]);
        subtitle(subtitle_str);
    
        xlabel(texlabel('D_{h} Multiple'));
        ylabel(texlabel('D_{v} Multiple'));
    
        hold on;
        plot(D_h_Multiples, ones(size(D_h_Multiples)), 'k-');
        plot(ones(size(D_h_Multiples)), D_v_Multiples, 'k-');
        hold off;
        
        out_path       = params.out_path;
        exp_path       = fullfile(out_path, params.exp_name);
    
        if ~(isfolder(exp_path))
            mkdir(exp_path);
        end
    
        addpath(out_path);
        addpath(exp_path);
    
        set(gcf, 'PaperUnits', 'inches');
        figWidth  = 7.5; % Figure width in inches.
        figHeight = 6; % Figure height in inches.
        set(gcf,...
            'PaperPosition', [0, 0, figWidth, figHeight],...
            'PaperSize', [figWidth, figHeight])
    
        file_name = sprintf('growthrates_Dv_Dh_%d_%d.pdf', ...
            wavenum, B_vs_idx);
        plot_file = fullfile(exp_path, file_name);
        print(plot_file, '-dpdf', '-painters');
    end
end