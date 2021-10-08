function lin_state = osh19_calc_lin_state(params, grid, bg_profs, mode, wavenum)

lin_state = struct();

% Extract some variables from structs to clean up later code
nx = params.nx;
ny = params.ny;
nz = params.nz;

lin_state.zeta_tau = zeros([ny, nx]);
lin_state.u_psi    = zeros([ny, nx, nz + 1]);
lin_state.v_psi    = zeros([ny, nx, nz + 1]);
lin_state.theta    = zeros([ny, nx, nz + 1]);
lin_state.q        = zeros([ny, nx, nz + 1]);
lin_state.tau_z    = zeros([ny, nx]);
lin_state.u_tau    = zeros([ny, nx]);
lin_state.v_tau    = zeros([ny, nx]);
lin_state.u        = zeros([ny, nx, nz + 1]);
lin_state.v        = zeros([ny, nx, nz + 1]);
lin_state.w        = zeros([ny, nx, nz + 1]);
lin_state.p        = zeros([ny, nx, nz + 1]);

H         = params.H;
beta      = params.beta;
g         = params.g;
c_p       = params.c_p;
L_v       = params.L_v;
B         = params.B;
theta_0   = params.theta_0;
tau_u     = params.tau_u;
tau_theta = params.tau_theta;
tau_up    = params.tau_up;
tau_mid   = params.tau_mid;

dz = grid.dz;

zz             = bg_profs.zz;
q_bg_vec_int   = bg_profs.q_bg_vec_int;
ddz_q_bg_vec   = bg_profs.ddz_q_bg_vec;
q_bg_merid_vec = bg_profs.q_bg_merid_vec;
tau_vec        = bg_profs.tau_vec;
D_h_vec        = bg_profs.D_h_vec;
D_v_vec        = bg_profs.D_v_vec;

k_x_dim = 2 * pi * wavenum / (2 * pi * params.P_E); % Dimensional zonal wavenumber

% Make Y matrix that you'll need
Y = diag(grid.yy);

% Make F (called Fmat to avoid confusion with flux vector F in nonlinear
% model)
F_mat     = fft(eye(ny))  / sqrt(ny);
F_mat_inv = ifft(eye(ny)) * sqrt(ny);

% Make spectral differentiation operators
spec_ddx   = 1i * k_x_dim * eye(ny);
spec_d2dx2 = -k_x_dim^2 * eye(ny);

k_y_vec     = [0:ny/2 -ny/2+1:-1];
k_y_vec_dim = 2 * pi * k_y_vec / ( 2 * params.P_Y);
spec_ddy    = diag(1i*k_y_vec_dim);
spec_ddy(ny/2+1,ny/2+1) = 0;
spec_d2dy2  = diag(-k_y_vec_dim.^2);
ddy         = F_mat_inv * spec_ddy * F_mat;
ddy_beta_y  = diag(ddy * beta * (grid.yy)');

spec_inv_lapl = inv(spec_d2dx2 + spec_d2dy2);


% Create block matrix A to represent linear system
neqns = (4 * nz - 3) * ny;
A=zeros(neqns,neqns);

% Populate A except for last two row and column blocks
for ii = 1:nz-1
    for jj = 1:nz-1

        % Populate the U parts of the matrix Aij
        Uuij = zeros(ny,ny);
        if ii == jj
            Uuij = Uuij + (1/tau_u) * eye(ny);
        end
        
        if ii == jj
            Uvij = -beta * F_mat * Y * F_mat_inv;
        else
            Uvij = zeros(ny,ny);
        end
        
        if ii > jj
            Utij =  g * dz * 1i * k_x_dim * jj        ...
                / (theta_0 * nz) * eye(ny);
        else
            Utij = -g * dz * 1i * k_x_dim * (nz - jj) ...
                / (theta_0 * nz) * eye(ny);
        end
        
        Uqij = zeros(ny,ny);
        
        % Populate the V parts of the matrix Aij
        if ii == jj
            Vuij = beta * F_mat * Y * F_mat_inv;
        else
            Vuij = zeros(ny,ny);
        end
        
        Vvij = zeros(ny,ny);
        if ii == jj
            Vvij = Vvij + (1/tau_u) * eye(ny);
        end
        
        if ii > jj
            Vtij =  g * dz * jj        ...
                / (theta_0 * nz) * spec_ddy;
        else
            Vtij = -g * dz * (nz - jj) ...
                / (theta_0 * nz) * spec_ddy;
        end
        
        Vqij = zeros(ny,ny);
        
        % Populate the Theta parts of the matrix Aij
        if ii >= jj
            Tuij = -1i * k_x_dim * B * dz * eye(ny);
        else
            Tuij = zeros(ny,ny);
        end
        
        if ii >= jj
            Tvij = -B * dz * spec_ddy;
        else
            Tvij = zeros(ny,ny);
        end
        
        Ttij = zeros(ny,ny);
        if ii == jj
            Ttij = Ttij + (1/tau_theta) * eye(ny);
        end
        
        mode2coeff = 0.5;
        if params.clin_conv_adj == 1
            Tqij = -sqrt(2) * L_v / c_p ...
                * ( 1/tau_up * ( sqrt(3) / nz * sin(zz(jj) * pi / H) ...
                                 - 4 * sqrt(3) / nz * sin( 2 * zz(jj) * pi / H) )...
                             * ( sin(zz(ii) * pi / H) ...
                                 - mode2coeff * sin(2 * zz(ii) * pi / H) ) ...
                    + 1/tau_mid * ( sqrt(3) / nz * sin(zz(jj) * pi / H) ...
                                    + 4 * sqrt(3) / nz * sin(2 * zz(jj) * pi / H) ) ...
                                * ( sin(zz(ii) * pi / H) ...
                                    + mode2coeff * sin(2 * zz(ii) * pi / H) ) ) ...
                * eye(ny);
        elseif params.clin_conv_adj == 2
            Tqij = -sqrt(2) * L_v / c_p ...
                * ( 1/tau_up * ( sqrt(3) / nz * sin(zz(jj) * pi / H) ...
                                 - 4 * sqrt(3) / nz * sin(2 * zz(jj) * pi / H) ) ...
                             * ( sin(zz(ii) * pi / H) ...
                                 - mode2coeff * sin(2 * zz(ii) * pi / H) ) ...
                    + 1/tau_mid * ( 2 / nz * sin(zz(jj) * pi / H) ) ...
                                * ( sin(zz(ii) * pi / H ) ) ) ...
                * eye(ny);
        else
            Tqij = zeros(ny,ny);
            if ii == jj
                Tqij = Tqij - L_v / c_p / tau_vec(ii) * eye(ny);
            end
        end
        
        
        % Populate the Q parts of the matrix Aij
        if ii >= jj
            Quij = -1i * k_x_dim * ddz_q_bg_vec(ii) * F_mat ...
                * diag(q_bg_merid_vec) * F_mat_inv * dz * eye(ny);
        else  
            Quij = zeros(ny,ny);
        end
        
        if ii >= jj
            Qvij = -ddz_q_bg_vec(ii) * F_mat * diag(q_bg_merid_vec) ...
                * F_mat_inv * dz * spec_ddy;
        else
            Qvij = zeros(ny,ny);
        end
        if ii == jj
            Qvij = Qvij + 1/2 * F_mat ...
                          * diag(ddy * q_bg_vec_int(ii) * transpose(q_bg_merid_vec)) ...
                          * F_mat_inv;
        end
        if ii == jj-1
            Qvij = Qvij + 1/2 * F_mat ...
                          * diag(ddy * q_bg_vec_int(ii) * transpose(q_bg_merid_vec)) ...
                          * F_mat_inv;
        end
        if ii == nz-1
            Qvij = Qvij - 1/2 * F_mat ...
                          * diag(ddy * q_bg_vec_int(ii) * transpose(q_bg_merid_vec)) ...
                          * F_mat_inv;
        end
        
        Qtij = zeros(ny,ny);
        if params.clin_conv_adj == 1
            Qqij = sqrt(2) ...
                * ( 1/tau_up * ( sqrt(3) / nz * sin(zz(jj) * pi / H) ...
                                 - 4 * sqrt(3) / nz * sin(2 * zz(jj) * pi / H) ) ...
                             * ( sin(zz(ii) * pi / H) ...
                                 - mode2coeff * sin(2 * zz(ii) * pi / H) ) ...
                    + 1/tau_mid * ( sqrt(3) / nz * sin(zz(jj) * pi / H) ...
                                    + 4 * sqrt(3) / nz * sin(2 * zz(jj) * pi / H) ) ...
                                * ( sin(zz(ii) * pi / H) ...
                                    + mode2coeff * sin(2 * zz(ii) * pi / H) ) ) ...
                * eye(ny);
            if ii == jj
                Qqij = Qqij + k_x_dim^2 * D_h_vec(ii) * eye(ny) ...
                    - D_h_vec(ii) * spec_d2dy2 ...
                    + 2 * D_v_vec(ii) / dz^2 * eye(ny);
            elseif ii == jj-1 || ii == jj+1
                Qqij = Qqij - D_v_vec(ii) / dz^2 * eye(ny);
            end
        elseif params.clin_conv_adj == 2
            Qqij = sqrt(2) ...
                * ( 1/tau_up * ( sqrt(3) / nz * sin(zz(jj) * pi / H) ...
                                 - 4 * sqrt(3) / nz * sin(2 * zz(jj) * pi / H) ) ...
                             * ( sin(zz(ii) * pi / H) ...
                                 - mode2coeff * sin(2 * zz(ii) * pi / H) ) ...
                    + 1/tau_mid * ( 2 / nz * sin(zz(jj) * pi / H) ) ...
                                * ( sin(zz(ii) * pi / H) ) ) ...
                * eye(ny);
            if ii == jj
                Qqij = Qqij + k_x_dim^2 * D_h_vec(ii) * eye(ny) ...
                    - D_h_vec(ii) * spec_d2dy2 ...
                    + 2 * D_v_vec(ii) / dz^2 * eye(ny);
            elseif ii == jj-1 || ii == jj+1
                Qqij = Qqij - D_v_vec(ii) / dz^2 * eye(ny);
            end
        else
            if ii == jj
                Qqij = 1/tau_vec(ii) * eye(ny) ...
                    + k_x_dim^2 * D_h_vec(ii) * eye(ny) ...
                    - D_h_vec(ii) * spec_d2dy2 ...
                    + 2*D_v_vec(ii) / dz^2 * eye(ny);
            elseif ii == jj-1 || ii == jj+1
                Qqij = -D_v_vec(ii) / dz^2 * eye(ny);
            else
                Qqij = zeros(ny,ny);
            end 
        end
        
        % Put all the parts together to make Aij, then populate A
        % All block rows and block columns except last two
        Aij = [Uuij Uvij Utij Uqij;...
               Vuij Vvij Vtij Vqij;...
               Tuij Tvij Ttij Tqij;...
               Quij Qvij Qtij Qqij];
        A(1 + (ii-1) * 4 * ny:ii * 4 * ny, ...
            1 + (jj-1) * 4 * ny:jj * 4 * ny) = Aij;
    end
end

% Last block row (not including last block column)
A(1 + (nz - 1) * 4 * ny:end, ...
    1:(nz - 1) * 4 * ny) = zeros(ny, 4 * ny * (nz - 1));

% Last block column (not including last block row)
for ii = 1:nz-1
    Uzij = zeros(ny,ny); 
    Vzij = zeros(ny,ny); 
    Tzij = zeros(ny,ny);
    Qzij = zeros(ny,ny);
    if ii == nz-1
        Qzij = Qzij + ...
            1/2 * F_mat * diag(ddy * q_bg_vec_int(ii) * transpose(q_bg_merid_vec)) ...
                * nz * F_mat_inv * spec_ddx * spec_inv_lapl;
    end
    Aij=[Uzij;...
         Vzij;...
         Tzij;...
         Qzij];
    A(1 + (ii - 1) * 4 * ny:ii * 4 * ny, 1 + (nz - 1) * 4 * ny:end) = Aij;
end

% Last block row and column
% Calculate derivative spectrally
Zzij = 1i * k_x_dim * F_mat * ddy_beta_y * F_mat_inv * spec_inv_lapl ...
    + 1/tau_u * eye(ny);

A(1 + (nz-1) * 4 * ny:end, 1 + (nz - 1) * 4 * ny:end) = Zzij;

% Find eigenvalues of A
[evecsFsp, evalsmat] = eig(A);
days_to_secs = 3600*24;
evals = -1i * diag(evalsmat) * days_to_secs; % evals now omega in 1/day units

sortbyspeed=0; 
% Positive - eastward;  negative - westward

% % Sort the eigenvalues according to speed/least damped. 
if sortbyspeed == 1 
    [~, Index_Sorted] = sort(real(evals),'descend');... % Fastest eastward speed first
else
    [~, Index_Sorted] = sort(imag(evals),'descend'); % largest growth rates first
end

evalssort = evals(Index_Sorted);

Growthrates = imag(evalssort);
disp(strcat(['Calculating linear mode #',int2str(mode),', k=',int2str(wavenumber)]));
disp('Growthrates');
disp(num2str(Growthrates(1:5)));

 
% Set sorted frequencies
Frequencies = real(evalssort) / (2*pi);

% Set sorted Right Evecs
evecsFspsort = evecsFsp(:, Index_Sorted);

% Now make a copy with the top-level winds replacing the barotropic
% relative vorticity
evecsphspYwinds = zeros(ny * (4 * nz - 3), 1);

% Inverse FFT to get back to physical space in y
for ii = 1:4*nz-3
    evecsphspYwinds(1 + (ii - 1) * ny:ii * ny, 1) = ...
        (F_mat_inv * evecsFspsort(1 + (ii-1) * ny:ii * ny, mode));
end

% Convert eigenvector to physical space in x
evectoshow = zeros(ny * (4 * nz - 3), nx);
for rowcount = 1:ny*(4*nz-3)
    evectoshow(rowcount,:) = real( evecsphspYwinds(rowcount,1) ...
                                   * exp(1i * k_x_dim * grid.xx) );
end

% Construct linear solution in physical space, diagnostic variables only
lin_state.zeta_tau = evectoshow(ny * 4 * (nz - 1) + 1:ny * 4 * (nz - 1) + ny, :);
for z_idx = 2:nz
    lin_state.u_psi(:,:,z_idx) = evectoshow(          1 + 4 * ny * (z_idx - 2):    ny + 4 * ny * (z_idx - 2), :);
    lin_state.v_psi(:,:,z_idx) = evectoshow(     ny + 1 + 4 * ny * (z_idx - 2):2 * ny + 4 * ny * (z_idx - 2), :);
    lin_state.theta(:,:,z_idx)  = evectoshow( 2 * ny + 1 + 4 * ny * (z_idx - 2):3 * ny + 4 * ny * (z_idx - 2), :);
    lin_state.q(:,:,z_idx)      = evectoshow( 3 * ny + 1 + 4 * ny * (z_idx - 2):4 * ny + 4 * ny * (z_idx - 2), :);
end

end
