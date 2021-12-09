function params_out = ens_convert_params(params_in)


days_to_secs = 3600*24;

params_out = params_in;

params_out.comm_freq = params_out.comm_freq * days_to_secs;

end

