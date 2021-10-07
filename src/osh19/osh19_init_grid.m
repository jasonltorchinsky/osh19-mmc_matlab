function grid = osh19_init_grid(params)

grid = struct();

% Set spatial grids
planet_eq_circ = 2 * pi * params.P_E;

dx = planet_eq_circ / params.nx;
dy = 2 * params.P_Y / params.ny;
dz = params.H / params.nz;

xx  = 0:dx:planet_eq_circ-dx;
yy  = -params.P_Y+dy/2:dy:params.P_Y-dy/2;
zzU = -dz/2:dz:params.H;
zzW = 0:dz:params.H;

[XX, YY]           = meshgrid(xx, yy);
[XX3U, YY3U, ZZ3U] = meshgrid(xx, yy, zzU);
[XX3W, YY3W, ZZ3W] = meshgrid(xx, yy, zzW);

grid.dx = dx;
grid.dy = dy;
grid.dz = dz;

grid.xx  = xx;
grid.yy  = yy;
grid.zzU = zzU;
grid.zzW = zzW;

grid.XX = XX;
grid.YY = YY;

grid.XX3U = XX3U;
grid.YY3U = YY3U;
grid.ZZ3U = ZZ3U;

grid.XX3W = XX3W;
grid.YY3W = YY3W;
grid.ZZ3W = ZZ3W;

% Set time-step size
c         = sqrt(params.g * params.B / params.theta_0) ...
    * params.H / pi;                             % Wave speed of the system (km s^(-1))
alpha_bar = (params.H / pi) * params.B;          % (K)
Q         = params.c_p * alpha_bar / params.L_v; % (Unitless)
F_scale   = (params.H / pi) * (params.B ...
    + params.L_v / params.c_p * params.B_vs ...
    - params.theta_0 * params.B_vs);             % (K)
F_tilde   = F_scale / alpha_bar;                 % (Unitless)                 
c_moist   = c * sqrt(F_tilde);                   % (km s^(-1))

CFL_cnst  = 0.15/40;                             % Not necessary?  
dt        = CFL_cnst * dx / c_moist;             % Initial dt, to be rounded
                                                 % to nearest fraction of day (s)

days_to_secs = 3600*24;
dt = dt / days_to_secs;
dt = ceil(1/dt);
dt = 1/dt * days_to_secs;

grid.dt = dt;

end

