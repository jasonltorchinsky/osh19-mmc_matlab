function state = cmg14_init_state(params)

state = struct();

if params.IC_type == 1
   % Zero initial condition
   
   state.u_1 = 3;
   state.u_2 = 3;
   state.v   = 0;
   state.w_u = 0;
    
end
            
end

