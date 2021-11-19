function [dcm_state_out, mjoo_state_out] = ens_comm(dcm_params, dcm_grid, ...
    dcm_state, mjoo_state, ens_params, eofs)

% Unpack common parameters
nx = dcm_params.nx;
nz = dcm_params.nz;

H_dcm = ens_params.dcm_comm.H;
B_dcm = ens_params.dcm_comm.B;
Lambda_dcm = ens_params.dcm_comm.Lambda;

H_mjoo = ens_params.mjoo_comm.H;
B_mjoo = ens_params.mjoo_comm.B;
Lambda_mjoo = ens_params.mjoo_comm.Lambda;

scale_x = (2 * pi) / (2 * pi * dcm_params.P_E);
scale_y = (2 * pi) / (2 * dcm_params.P_Y);

% Project DCM state onto non-dimensional concatenated EOF to get DCM MJO indices
u_proj_pri = ens_proj_u(dcm_params, dcm_grid, dcm_state, eofs);
q_proj_pri = ens_proj_q(dcm_params, dcm_grid, dcm_state, eofs);
a_proj_pri = [u_proj_pri q_proj_pri];

% Get MJO indices, climate parameter from DCM
u1_dcm_pri = a_proj_pri * eofs.raw_eof1;
u2_dcm_pri = a_proj_pri * eofs.raw_eof2;
v_dcm_pri  = dcm_params.B_vs;

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

mjoo_ptruth = [u1_mjoo_pri; u2_mjoo_pri; v_mjoo_pri];
mjoo_stoch = normrnd(0, 1, [1, 1]);
mjoo_obs = H_mjoo * mjoo_ptruth ...
    + Lambda_mjoo * mjoo_stoch;

% Communicate the MJO indices, climate parameter
gain_dcm = B_dcm * transpose(H_dcm) / (Lambda_dcm + H_dcm * B_dcm * transpose(H_dcm));
dcm_pst = dcm_pri - gain_dcm * (H_dcm * dcm_pri - dcm_obs);
% TRUST MJOO MODEL ENTIRELY
dcm_pst = mjoo_pri;

gain_mjoo = B_mjoo * transpose(H_mjoo) / (Lambda_mjoo + H_mjoo * B_mjoo * transpose(H_mjoo));
mjoo_pst = mjoo_pri - gain_mjoo * (H_mjoo * mjoo_pri - mjoo_obs);

% Update DCM state
% Get posterior projections
u1_dcm_pst = dcm_pst(1);
u2_dcm_pst = dcm_pst(2);
a_proj_pst = a_proj_pri + (u1_dcm_pst - u1_dcm_pri) * eofs.raw_eof1 ...
    + (u2_dcm_pst - u2_dcm_pri) * eofs.raw_eof2;
u_proj_pst = a_proj_pst(1:nx);
q_proj_pst = a_proj_pst(nx+1:end);

% Un-project the prior and posterior u, q, and update them in the state
u_pri = ens_unproj_u(dcm_params, dcm_grid, u_proj_pri, eofs);
q_pri = ens_unproj_q(dcm_params, dcm_grid, q_proj_pri, eofs);

u_pst = ens_unproj_u(dcm_params, dcm_grid, u_proj_pst, eofs);
q_pst = ens_unproj_q(dcm_params, dcm_grid, q_proj_pst, eofs);

dcm_state_out = dcm_state;
dcm_state_out.q = dcm_state_out.q - q_pri + q_pst;
% u is a prognostic variable, so we need to correct the diagnostic variables too
dcm_state_out.u = dcm_state_out.u - u_pri + u_pst;
dcm_state_out.u(:,:,1)    = 0.0; % Enforce BCs
dcm_state_out.u(:,:,nz+1) = 0.0;
dcm_state_out.u_tau = squeeze(mean(dcm_state_out.u, 3));
for kk = 1:nz+1
    dcm_state_out.u_psi(:,:,kk) = dcm_state_out.u(:,:,kk) - dcm_state_out.u_tau;
end
dcm_state_out.zeta_tau = D1(dcm_state_out.v_tau, 'x', scale_x) ...
    - D1(dcm_state_out.u_tau, 'y', scale_y);
dcm_state_out = osh19_prognose_state(dcm_params, dcm_grid, dcm_state_out);

% Update MJOO state
mjoo_state_out = mjoo_state;
%mjoo_state_out.u_1 = mjoo_pst(1);
%mjoo_state_out.u_2 = mjoo_pst(2);
%mjoo_state_out.v   = mjoo_pst(3);

fprintf('u1_dcm_pri,  u2_dcm_pri:  %.4f, %.4f\n',   u1_dcm_pri,  u2_dcm_pri);
fprintf('u1_mjoo_pri, u2_mjoo_pri: %.4f, %.4f\n\n', u1_mjoo_pri, u2_mjoo_pri);

fprintf('u1_dcm_pst,  u2_dcm_pst:  %.4f, %.4f\n',   u1_dcm_pst,  u2_dcm_pst);
fprintf('u1_mjoo_pst, u2_mjoo_pst: %.4f, %.4f\n\n', mjoo_pst(1), mjoo_pst(2));

fprintf('Max u_proj in, max u_proj out: %.4f, %.4f\n\n', ...
    max(abs(u_pri), [], 'all'), max(abs(u_pst), [], 'all'));

fprintf('Max u in, max u out: %.4f, %.4f\n\n', ...
    max(abs(dcm_state.u), [], 'all'), max(abs(dcm_state_out.u), [], 'all'));
end

