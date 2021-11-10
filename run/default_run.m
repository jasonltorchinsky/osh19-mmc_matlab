out_path = 'output';
exp_name = 'long_default';

sim_days = 400;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters for the truth (OSH19) [Default]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
truth_params.sim_days = sim_days; % Number of days to simulate (d)
truth_params.out_freq = 1; % How often to output data (d)

truth_params.out_path = out_path;
truth_params.exp_name = exp_name;
truth_params.component_name = 'truth';

truth_params.init_simulation = false;
truth_params.run_simulation = false;
truth_params.create_plots = false;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters for the Empirical Orthogonal Function calculation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eof_params = struct();

eof_params.Q_mode = 'up'; % Use Q_mid ('mid') or Q_up ('up')

eof_params.out_path = out_path;
eof_params.exp_name = exp_name;
eof_params.component_name = truth_params.component_name;

eof_params.calc_eofs = false;
eof_params.create_plots = false;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters for the deficient climate model (OSH19) [Case 2]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dcm_params = struct();

dcm_params.nx = 32; % Number of points zonally
dcm_params.ny = 32; % Number of points meridionally
dcm_params.nz = 8;  % Number of non-trivial levels for u, v, p
dcm_params.H = 16;  % Height of troposphere (km)
dcm_params.P_Y = 6000; % Distance from equator to channel wall (km)
dcm_params.P_E = 40000/(2*pi); % Radius of planet (km)
dcm_params.beta = 2.3 * 10^(-8); % Meridional variation of Coriolis parameter (km^(-1) s^(-1))
dcm_params.g = 9.8 * 10^(-3); % Acceleration due to gravity (km s^(-2))
dcm_params.c_p = 1000; % Specific heat of dry air at constant pressure (J kg^(-1) K^(-1))
dcm_params.L_v = 2.5 * 10^(6); % Latent heat of vaporization (J kg^(-1))
dcm_params.B = 3; % Background potential temperature vertical gradient (K km^(-1))
dcm_params.theta_0 = 300; % Reference potential temperature (K)
dcm_params.tau_u = 16; % Wind damping time-scale (d)
dcm_params.tau_theta = 16; % Potential temperature damping time-scale (d)
dcm_params.tau_up = 1; % Moisture damping time-scale in upper-troposphere (d)
dcm_params.tau_mid = 4/24; % Moisture damping time-scale in mid-troposphere (d)
dcm_params.B_vs = -1.28 * 10^(-3); % Mean vertical q_bg gradient (kg kg^(-1) km^(-1))
dcm_params.a = 0.5; % 1--pole-to-equator q_bg ratio
dcm_params.L_tilde = 2000; % q_bg meridional decay length scale (km)
dcm_params.s_tilde = 12; % q_bg vertical decay length scale (km)
dcm_params.D_hUp = 121.6; % Horizontal q diffusion in upper-troposphere (km^(2) s^(-1))
dcm_params.D_hMid = 15.2; % Horizontal q diffusion in mid-troposphere (km^(2) s^(-1))
dcm_params.D_v = 5 * 10^(-4); % Vertical q diffusion (km^(2) s^(-1))
dcm_params.IC_type = 1; % Initial condition type:
                          % 1 = get modes, wavenumbers from linear solution
                          % 2 = load state from file
dcm_params.IC_modes = [1]; % Modes to use for initial condition (IC)
dcm_params.IC_wavenums = [1, 2, 3]; % Zonal wavenumber to use for IC
dcm_params.IC_amp = 1; % Amplification factor for IC
dcm_params.clin_conv_adj = 1; % Options for baroclinic modes in IC
dcm_params.sim_days = sim_days; % Number of days to simulate (d)
dcm_params.out_freq = 1; % How often to output data (d)

dcm_params.out_path = out_path;
dcm_params.exp_name = exp_name;
dcm_params.component_name = 'dcm';

dcm_params.init_simulation = false;
dcm_params.run_simulation = false;
dcm_params.create_plots = false;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters for the MJO-only model (OSH19) [Default]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mjoo_params = struct();

mjoo_params.d_u     = 0.5*0.9; % Damping for MJO modes (m^(-1))
mjoo_params.d_v     = 0.9; % Damping for stochastic damping (m^(-1))
mjoo_params.d_w     = 0.5; % Damping for stochastic phase (m^(-1))
mjoo_params.gamma   = 0.3; % Strength of non-linear interaction (m^(-1))
mjoo_params.a       = 1.5; % Background state phase of MJO modes (m^(-1))
mjoo_params.w_u_hat = 0; % Background mean state of stochastic phase (m^(-1))

mjoo_params.sigma_u = 0.3; % Strength of stochastic forcing for MJO modes (m^(-1/2))
mjoo_params.sigma_v = 1; % Strength of stochastic forcing for stochastic damping (m^(-1/2))
mjoo_params.sigma_w = 1.1; % Strength of stochastic forcing for stochastic phase (m^(-1/2))

mjoo_params.f_0     = 1; % Mean time-periodic damping
mjoo_params.f_t     = 6.9; % Amplitude of time-periodic damping
mjoo_params.w_f     = 2 * pi / 12; % Frequency of time-periodic damping (m^(-1))
mjoo_params.phi     = -1; % Phase-shift of time-periodic damping

mjoo_params.dt      = 1/30; % Time-step size (m) !! WATCH THIS WITH OUTPUT FREQUENCY !!
                            %                    !! MUST BE LESS THAN OUTPUT
                            %                       FREQUENCY !!

mjoo_params.IC_type = 1; % Initial condition type:
                         % 1 = zero
                         % 2 = load state from file

mjoo_params.sim_days = sim_days; % Number of days to simulate (d)
mjoo_params.out_freq = 1; % How often to output data (d)

mjoo_params.out_path = out_path;
mjoo_params.exp_name = exp_name;
mjoo_params.component_name = 'mjoo';

mjoo_params.init_simulation = true;
mjoo_params.run_simulation = true;
mjoo_params.create_plots = true;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters for the Multi-Model Ensemble
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ens_params = struct();

% Ensemble parameters
ens_params.out_path = out_path;
ens_params.exp_name = exp_name;
ens_params.component_name = 'ens';

ens_params.init_simulation = false;
ens_params.run_simulation = false;
ens_params.create_plots = false;

% Ensemble DCM parameters
ens_params.dcm_params = dcm_params;

ens_params.dcm_params.component_name = 'ens_dcm';

ens_params.dcm_params.init_simulation = ens_params.init_simulation;
ens_params.dcm_params.run_simulation = ens_params.run_simulation;
ens_params.dcm_params.create_plots = ens_params.create_plots;

% Ensemble MJOO parameters
ens_params.mjoo_params = mjoo_params;

ens_params.mjoo_params.component_name = 'ens_mjoo';

ens_params.mjoo_params.init_simulation = ens_params.init_simulation;
ens_params.mjoo_params.run_simulation = ens_params.run_simulation;
ens_params.mjoo_params.create_plots = ens_params.create_plots;

% Ensemble EOF parameters
ens_params.eof_params = eof_params;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Miscellaneous Parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

misc_params = struct();

misc_params.create_plots = true;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Run the simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

error = main(truth_params, eof_params, dcm_params, mjoo_params, ens_params);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot truth simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot EOFs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if eof_params.create_plots
    clf('reset');
    
    eof_plot_eofs(eof_params, 10.0);
    
    clf('reset');
    
    eof_plot_exp(eof_params);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot DCM simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if dcm_params.create_plots
    clf('reset');
    
    osh19_plot_growths_freqs(dcm_params, 1);
    
    clf('reset');
    
    osh19_plot_evo(dcm_params, 10.0, 0.0);
    
    clf('reset');
    
    osh19_plot_hovmoller(dcm_params, 10.0, 0.0);
    
    if dcm_params.sim_days > 256
        clf('reset');
    
        osh19_plot_wheeler_kiladis(dcm_params, 128, 38, 0.1, 12);
    end
    
    clf('reset');
    
    osh19_plot_energy_evo(dcm_params, 12);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot MJOO simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if mjoo_params.create_plots
    clf('reset');
    
    cmg14_plot_evo(mjoo_params);
    
    clf('reset');
    
    cmg14_plot_pdfs(mjoo_params);
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot ensemble simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ens_params.create_plots
    clf('reset');
    
    osh19_plot_growths_freqs(ens_params.dcm_params, 1);
    
    clf('reset');
    
    osh19_plot_evo(ens_params.dcm_params, 10.0, 0.0);
    
    clf('reset');
    
    osh19_plot_hovmoller(ens_params.dcm_params, 10.0, 0.0);
    
    if ens_params.dcm_params.sim_days > 256
        clf('reset');
    
        osh19_plot_wheeler_kiladis(ens_params.dcm_params, 128, 38, 0.1, 12);
    end
    
    clf('reset');
    
    osh19_plot_energy_evo(ens_params.dcm_params, 12);
    
    clf('reset');
    
    cmg14_plot_evo(ens_params.mjoo_params);
    
    clf('reset');
    
    cmg14_plot_pdfs(ens_params.mjoo_params);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot miscellaneous figures
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if misc_params.create_plots
   clf('reset');
   
   misc_plot_mjos(eof_params, mjoo_params, 10.0, 0.0)
end