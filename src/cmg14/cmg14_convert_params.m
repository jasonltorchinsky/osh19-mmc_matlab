function params_out = cmg14_convert_params(params_in)

months_to_days = 30;
days_to_secs   = 3600 * 24;
months_to_secs = months_to_days * days_to_secs;

params_out = params_in;

params_out.d_u = params_out.d_u / months_to_secs;
params_out.d_v = params_out.d_v / months_to_secs;
params_out.d_w = params_out.d_w / months_to_secs;

params_out.a       = params_out.a / months_to_secs;
params_out.w_u_hat = params_out.w_u_hat / months_to_secs;

params_out.sigma_u = params_out.sigma_u / sqrt(months_to_secs);
params_out.sigma_v = params_out.sigma_v / sqrt(months_to_secs);
params_out.sigma_w = params_out.sigma_w / sqrt(months_to_secs);

params_out.f_0 = params_out.f_0 / months_to_secs;
params_out.f_t = params_out.f_t / months_to_secs;
params_out.w_f = params_out.w_f / months_to_secs;

params_out.dt  = params_out.dt * months_to_secs;

params_out.sim_days = params_out.sim_days * days_to_secs;
params_out.out_freq = params_out.out_freq * days_to_secs;
end

