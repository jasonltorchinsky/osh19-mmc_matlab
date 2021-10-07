function ierror = main(osh19_params)

ierror = 1;

%osh19_grid     = struct();
%osh19_bg_profs = struct();
%osh19_state    = struct();

osh19_params   = osh19_convert_params(osh19_params);
osh19_grid     = osh19_init_grid(osh19_params);
osh19_bg_profs = osh19_init_bg_profs(osh19_params, osh19_grid);
osh19_state    = osh19_init_state(osh19_params, osh19_grid, osh19_bg_profs);

disp(osh19_state);

end