function eofs = eof_calc_eof(params)

eofs = struct();

% Get vertical, meridional projections of u, q 
u_proj = eof_proj_u(params);
q_proj = eof_proj_q(params);

% Take out time-series means
u_proj = detrend(u_proj, 0);
q_proj = detrend(q_proj, 0);

% Normalize u_proj, q_proj by their variance
u_proj_var = var(u_proj, 0, 'all');
u_proj = u_proj / u_proj_var;

q_proj_var = var(q_proj, 0, 'all');
q_proj = q_proj / q_proj_var;

% Perform the EOF decompostion via eigenvalue decomposition
comb_proj = [u_proj q_proj]; % Combined projection matrix
cov = transpose(comb_proj) * comb_proj;
[eof_mtx, L] = eig(cov);

% Get the first two EOFs, and expansion coefficients
[nt, nx] = size(u_proj);

eofs.u_eof1 = transpose(eof_mtx(1:nx, end)) * u_proj_var;
eofs.u_eof2 = transpose(eof_mtx(1:nx, end-1)) * u_proj_var;
eofs.q_eof1 = transpose(eof_mtx(nx+1:end, end)) * q_proj_var;
eofs.q_eof2 = transpose(eof_mtx(nx+1:end, end-1)) * q_proj_var;

exp_coeffs = comb_proj * eof_mtx(:, [end end-1]);
eofs.exp1 = transpose(exp_coeffs(1:nt, end));
eofs.exp2 = transpose(exp_coeffs(1:nt, end-1));

eofs.scf = sort(diag(L) / trace(L), 'descend'); % Scaled covariance fraction

eofs = eof_recon_mjo(params, eofs);

end

