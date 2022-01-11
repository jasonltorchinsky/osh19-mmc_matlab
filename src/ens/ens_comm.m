function [dcm_state_out, mjoo_state_out] = ens_comm(dcm_params, dcm_grid, ...
    dcm_bg_profs, dcm_state, mjoo_state, ens_params, eofs)

% Unpack common parameters
nx = dcm_params.nx;
ny = dcm_params.ny;
nz = dcm_params.nz;

dz = dcm_grid.dz;

H_dcm = ens_params.dcm_comm.H;
B_dcm = ens_params.dcm_comm.B;
Lambda_dcm = ens_params.dcm_comm.Lambda;

H_mjoo = ens_params.mjoo_comm.H;
B_mjoo = ens_params.mjoo_comm.B;
Lambda_mjoo = ens_params.mjoo_comm.Lambda;

% Project DCM state onto non-dimensional concatenated EOF to get DCM MJO indices
u_proj_pri = ens_proj_u(dcm_params, dcm_grid, dcm_state.u, eofs);
q_proj_pri = ens_proj_q(dcm_params, dcm_grid, dcm_state.q, eofs);
a_proj_pri = [u_proj_pri q_proj_pri];

% Get mean vertical moisture gradient.
q_tot = dcm_state.q + dcm_bg_profs.q_bg_mat;
ddz_q_tot = zeros([nx, ny, nz-1]);
for kk = 2:nz
    ddz_q_tot(:,:,kk-1) = (q_tot(:,:,kk+1) - q_tot(:,:,kk-1)) / (2 * dz);
end
mean_ddz_q_tot = mean(ddz_q_tot, 'all');

% Get MJO indices, climate parameter from DCM
u1_dcm_pri = a_proj_pri * eofs.eof1;
u2_dcm_pri = a_proj_pri * eofs.eof2;
v_dcm_pri  = mean_ddz_q_tot/dcm_params.B_vs;

dcm_pri = [u1_dcm_pri; u2_dcm_pri; v_dcm_pri];

% Get MJO indices, climate parameter from MJOO
u1_mjoo_pri = mjoo_state.u_1;
u2_mjoo_pri = mjoo_state.u_2;
v_mjoo_pri  = mjoo_state.v;

mjoo_pri = [u1_mjoo_pri; u2_mjoo_pri; v_mjoo_pri];

% Generate "observations"
dcm_ptruth = [u1_mjoo_pri; u2_mjoo_pri; v_dcm_pri];
dcm_stoch = normrnd(0, 1, [2, 1]);
dcm_obs = H_dcm * dcm_ptruth + Lambda_dcm * dcm_stoch;

mjoo_ptruth = [u1_mjoo_pri; u2_mjoo_pri; v_dcm_pri];
mjoo_stoch = normrnd(0, 1, [1, 1]);
mjoo_obs = H_mjoo * mjoo_ptruth ...
    + Lambda_mjoo * mjoo_stoch;

% Communicate the MJO indices, climate parameter
gain_dcm = B_dcm * transpose(H_dcm) / (Lambda_dcm + H_dcm * B_dcm * transpose(H_dcm));
dcm_pst = dcm_pri - gain_dcm * (H_dcm * dcm_pri - dcm_obs);
% Manually set DCM parameters
%dcm_pst(1:2) = mjoo_pri(1:2);
%dcm_pst(3)   = dcm_pri(3);


gain_mjoo = B_mjoo * transpose(H_mjoo) / (Lambda_mjoo + H_mjoo * B_mjoo * transpose(H_mjoo));
mjoo_pst = mjoo_pri - gain_mjoo * (H_mjoo * mjoo_pri - mjoo_obs);
% Manually set MJOO parameters
%mjoo_pst(1:2) = mjoo_pri(1:2);
%mjoo_pst(3)   = dcm_pri(3);

% Update DCM state
u1_dcm_pst = dcm_pst(1);
u2_dcm_pst = dcm_pst(2);
delta_mjo_proj = (u1_dcm_pst - u1_dcm_pri) * eofs.eof1.' ...
    + (u2_dcm_pst - u2_dcm_pri) * eofs.eof2.';
delta_u_mjo_proj = delta_mjo_proj(1:nx);
delta_q_mjo_proj = delta_mjo_proj(nx+1:end);

% % Un-project the prior and posterior u, q, and update them in the state
delta_u_mjo = ens_unproj_u(dcm_params, dcm_grid, delta_u_mjo_proj, eofs);
delta_q_mjo = ens_unproj_q(dcm_params, dcm_grid, delta_q_mjo_proj, eofs);


dcm_state_out = dcm_state;
dcm_state_out.u = dcm_state_out.u + delta_u_mjo;
dcm_state_out.q = dcm_state_out.q + delta_q_mjo;

% u is a prognostic variable, so we need to correct the diagnostic variables too
scale_x = (2 * pi) / (2 * pi * dcm_params.P_E);
scale_y = (2 * pi) / (2 * dcm_params.P_Y);
% Update u_psi, w
for kk = 2:nz
    dcm_state_out.u_psi(:,:,kk) = dcm_state_out.u(:,:,kk) - dcm_state_out.u_tau;
end
for kk = 2:nz
    dcm_state_out.w(:,:,kk) = dcm_state_out.w(:,:,kk-1) ...
        - dz * (D1(dcm_state_out.u(:,:,kk),'x',scale_x) ...
        + D1(dcm_state_out.v(:,:,kk),'y',scale_y));
end

% Update MJOO state
mjoo_state_out = mjoo_state;
mjoo_state_out.u_1 = mjoo_pst(1);
mjoo_state_out.u_2 = mjoo_pst(2);
mjoo_state_out.v   = mjoo_pst(3);


% WARNING: BCS ARE NOT ENFORCED

% % Check if project and un-project are really inverses
% a_mjo_proj = 0 * eofs.eof1.' + 1 * eofs.eof2.';
% u_mjo_proj = a_mjo_proj(1:nx);
% q_mjo_proj = a_mjo_proj(nx+1:end);
% 
% u_mjo = ens_unproj_u(dcm_params, dcm_grid, u_mjo_proj, eofs);
% q_mjo = ens_unproj_q(dcm_params, dcm_grid, q_mjo_proj, eofs);
% 
% u_mjo_reproj = ens_proj_u(dcm_params, dcm_grid, u_mjo, eofs);
% [q_mjo_reproj, ~, ~] = ens_proj_q(dcm_params, dcm_grid, q_mjo, eofs);
% a_mjo_reproj = [u_mjo_reproj q_mjo_reproj];
% 
% u1 = a_mjo_reproj * eofs.eof1;
% u2 = a_mjo_reproj * eofs.eof2;

% % Take dcm_state_out.u and get the MJO indices from there!
% u_out_proj = ens_proj_u(dcm_params, dcm_grid, dcm_state_out.u, eofs);
% q_out_proj = ens_proj_q(dcm_params, dcm_grid, dcm_state_out.q, eofs);
% a_out_proj = [u_out_proj q_out_proj];
% 
% % Get MJO indices, climate parameter from DCM
% u1_out_dcm = a_out_proj * eofs.eof1;
% u2_out_dcm = a_out_proj * eofs.eof2;
% v_out_dcm  = dcm_params.B_vs;
% 
% disp([u1_out_dcm, u2_out_dcm]);

% % Check out the prior and posterior MJOs
% a_mjo_proj_pri = u1_dcm_pri * eofs.eof1.' + u2_dcm_pri * eofs.eof2.';
% u_mjo_proj_pri = a_mjo_proj_pri(1:nx);
% q_mjo_proj_pri = a_mjo_proj_pri(nx+1:end);
% 
% u_mjo_pri = ens_unproj_u(dcm_params, dcm_grid, u_mjo_proj_pri, eofs);
% q_mjo_pri = ens_unproj_q(dcm_params, dcm_grid, q_mjo_proj_pri, eofs);
% 
% a_mjo_proj_pst = u1_dcm_pst * eofs.eof1.' + u2_dcm_pst * eofs.eof2.';
% u_mjo_proj_pst = a_mjo_proj_pst(1:nx);
% q_mjo_proj_pst = a_mjo_proj_pst(nx+1:end);
% 
% u_mjo_pst = ens_unproj_u(dcm_params, dcm_grid, u_mjo_proj_pst, eofs);
% q_mjo_pst = ens_unproj_q(dcm_params, dcm_grid, q_mjo_proj_pst, eofs);

% fprintf('u1_dcm_pri,  u2_dcm_pri:  %.4f, %.4f\n',   u1_dcm_pri,  u2_dcm_pri);
% fprintf('u1_mjoo_pri, u2_mjoo_pri: %.4f, %.4f\n\n', u1_mjoo_pri, u2_mjoo_pri);
% 
% fprintf('u1_dcm_pst,  u2_dcm_pst:  %.4f, %.4f\n',   u1_dcm_pst,  u2_dcm_pst);
% fprintf('u1_mjoo_pst, u2_mjoo_pst: %.4f, %.4f\n\n', mjoo_pst(1), mjoo_pst(2));
% 
% fprintf('Max u_proj in, max u_proj out: %.4f, %.4f\n\n', ...
%     max(abs(u_pri), [], 'all'), max(abs(u_pst), [], 'all'));

% fprintf('Max u in, max u out: %.4f, %.4f\n\n', ...
%    max(abs(dcm_state.u), [], 'all')*10^3, ...
%    max(abs(dcm_state_out.u), [], 'all')*10^3);
% 
% x = max(abs(dcm_state.u - dcm_state_out.u), [], 'all') * 10^3;
% y = max(abs(dcm_state.u), [], 'all') * 10^3;
% fprintf('Max u change: %.8f \n\n', x);


end

