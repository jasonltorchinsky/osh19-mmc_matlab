out_path = 'output';
exp_name = 'default';

% Parameters for the truth (OSH19)
truth_params = struct();

truth_params.nx = 32; % Number of points zonally
truth_params.ny = 32; % Number of points meridionally
truth_params.nz = 8;  % Number of non-trivial levels for u, v, p
truth_params.H = 16;  % Height of troposphere (km)
truth_params.P_Y = 6000; % Distance from equator to channel wall (km)
truth_params.P_E = 40000/(2*pi); % Radius of planet (km)
truth_params.beta = 2.3 * 10^(-8); % Meridional variation of Coriolis parameter (km^(-1) s^(-1))
truth_params.g = 9.8 * 10^(-3); % Acceleration due to gravity (km s^(-2))
truth_params.c_p = 1000; % Specific heat of dry air at constant pressure (J kg^(-1) K^(-1))
truth_params.L_v = 2.5 * 10^(6); % Latent heat of vaporization (J kg^(-1))
truth_params.B = 3; % Background potential temperature vertical gradient (K km^(-1))
truth_params.theta_0 = 300; % Reference potential temperature (K)
truth_params.tau_u = 25; % Wind damping time-scale (d)
truth_params.tau_theta = 25; % Potential temperature damping time-scale (d)
truth_params.tau_up = 1; % Moisture damping time-scale in upper-troposphere (d)
truth_params.tau_mid = 2/24; % Moisture damping time-scale in mid-troposphere (d)
truth_params.B_vs = -1.34 * 10^(-3); % Mean vertical q_bg gradient (kg kg^(-1) km^(-1))
truth_params.a = 0.25; % 1--pole-to-equator q_bg ratio
truth_params.L_tilde = 2000; % q_bg meridional decay length scale (km)
truth_params.s_tilde = 12 * 10^(4); % q_bg vertical decay length scale (km)
truth_params.D_hUp = 60.8; % Horizontal q diffusion in upper-troposphere (km^(2) s^(-1))
truth_params.D_hMid = 7.6; % Horizontal q diffusion in mid-troposphere (km^(2) s^(-1))
truth_params.D_v = 1 * 10^(-4); % Vertical q diffusion (km^(2) s^(-1))
truth_params.IC_type = 1; % Initial condition type:
                          % 1 = get modes, wavenumbers from linear solution
                          % 2 = load state from file
truth_params.IC_modes = [1]; % Modes to use for initial condition (IC)
truth_params.IC_wavenums = [1, 2, 3]; % Zonal wavenumber to use for IC
truth_params.IC_amp = 1; % Amplification factor for IC
truth_params.clin_conv_adj = 2; % Options for baroclinic modes in IC
truth_params.sim_days = 100; % Number of days to simulate (d)
truth_params.out_freq = 1; % How often to output data (d)

truth_params.out_path = out_path;
truth_params.exp_name = exp_name;
truth_params.component_name = 'truth';

truth_params.init_simulation = false;
truth_params.run_simulation = false;
truth_params.create_plots = false;

% Parameters for the Empirical Orthogonal Function calculation
eof_params = struct();

eof_params.Q_mode = 'mid'; % Use Q_mid ('mid') or Q_up ('up')

eof_params.out_path = out_path;
eof_params.exp_name = exp_name;
eof_params.component_name = truth_params.component_name;

eof_params.calc_eofs = true;
eof_params.create_plots = false;

% Run the simulation and plotting

error = main(truth_params, eof_params);

if truth_params.create_plots
    clf('reset');
    
    osh19_plot_growths_freqs(truth_params, 1);
    
    clf('reset');
    
    osh19_plot_evo(truth_params, 10.0, 0.0);
    
    clf('reset');
    
    osh19_plot_hovmoller(truth_params, 10.0, 0.0);
    
    if truth_params.sim_days > 256
        clf('reset');
    
        osh19_plot_wheeler_kiladis(truth_params, 128, 38, 0.1, 12);
    end
    
    clf('reset');
    
    osh19_plot_energy_evo(truth_params, 12);
end

if eof_params.create_plots
    clf('reset');
    
    eof_plot_eofs(eof_params, 10.0);
    
    clf('reset');
    
    eof_plot_exp(eof_params);
end
    