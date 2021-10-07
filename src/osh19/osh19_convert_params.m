function params_out = osh19_convert_params(params_in)

days_to_secs = 3600*24;

params_out = params_in;

params_out.tau_u     = params_in.tau_u     * days_to_secs;
params_out.tau_theta = params_in.tau_theta * days_to_secs;
params_out.tau_up    = params_in.tau_up    * days_to_secs;
params_out.tau_mid   = params_in.tau_mid   * days_to_secs;

end

