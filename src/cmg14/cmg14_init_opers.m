function opers = cmg14_init_opers(params)

% unpack some paramaters to simplify equations
d_u = params.d_u;
d_v = params.d_v;
d_w = params.d_w;
a   = params.a;

dt  = params.dt;

% The linear-part operator is constant throughout the simulation, so we store it
A = [[-d_u -a 0 0]
    [a -d_u 0 0]
    [0 0 -d_v 0]
    [0 0 0 -d_w]];
opers.A = inv(eye(4) - dt*A);


% The non-linear-part, deterministic forcing-part, and stochastic forcing-part
% change throuhgout simulation, but we'll keep them here in case we need them
opers.B = zeros([4,1]); % Non-linear-part
opers.F = zeros([4,1]); % Determinstic forcing-part

end

