osh19_params = struct();
osh19_params.nx = 32; % Number of points zonally
osh19_params.ny = 32; % Number of points meridionally
osh19_params.nz = 8;  % Number of non-trivial levels for u, v, p
osh19_params.H = 16;  % Height of troposphere (km)
osh19_params.P_Y = 6000; % Distance from equator to channel wall (km)
osh19_params.P_E = 40000/(2*pi); % Radius of planet (km)
osh19_params.beta = 2.3 * 10^(-8); % Meridional variation of Coriolis parameter (km^(-1) s^(-1))
osh19_params.g = 9.8 * 10^(-3); % Acceleration due to gravity (km s^(-2))
osh19_params.c_p = 1000; % Specific heat of dry air at constant pressure (J kg^(-1) K^(-1))
osh19_params.L_v = 2.5 * 10^(6); % Latent heat of vaporization (J kg^(-1))
osh19_params.B = 3; % Background potential temperature vertical gradient (K km^(-1))
osh19_params.theta_0 = 300; % Reference potential temperature (K)
osh19_params.tau_u = 25; % Wind damping time-scale (d)
osh19_params.tau_theta = 25; % Potential temperature damping time-scale (d)
osh19_params.tau_up = 1; % Moisture damping time-scale in upper-troposphere (d)
osh19_params.tau_mid = 2/24; % Moisture damping time-scale in mid-troposphere (d)
osh19_params.B_vs = -1.34 * 10^(-3); % Mean vertical q_bg gradient (kg kg^(-1) km^(-1))
osh19_params.a = 0.25; % 1--pole-to-equator q_bg ratio
osh19_params.L_tilde = 2000; % q_bg meridional decay length scale (km)
osh19_params.s_tilde = 12 * 10^(4); % q_bg vertical decay length scale (km)
osh19_params.D_hUp = 60.8; % Horizontal q diffusion in upper-troposphere (km^(2) s^(-1))
osh19_params.D_hMid = 7.6; % Horizontal q diffusion in mid-troposphere (km^(2) s^(-1))
osh19_params.D_v = 1 * 10^(-4); % Vertical q diffusion (km^(2) s^(-1))
osh19_params.IC_type = 1; % Initial condition type:
                          % 1 = get modes, wavenumbers from linear solution
                          % 2 = load state from file
osh19_params.IC_modes = [1]; % Modes to use for initial condition (IC)
osh19_params.IC_wavenums = [1, 2, 3]; % Zonal wavenumber to use for IC
osh19_params.IC_amp = 1; % Amplification factor for IC
osh19_params.clin_conv_adj = 2; % Options for baroclinic modes in IC
osh19_params.sim_days = 100; % Number of days to simulate (d)
osh19_params.out_freq = 1; % How often to output data (d)

osh19_params.out_path = 'output';
osh19_params.exp_name = 'default';
osh19_params.component_name = 'truth';

osh19_params.run_simulation = true;
osh19_params.create_plots = true;

if osh19_params.run_simulation
    error = main(osh19_params);
end

if osh19_params.create_plots
    osh19_plot_evo(osh19_params, 9.0, 0.0)
end