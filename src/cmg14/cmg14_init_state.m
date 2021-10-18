function state = cmg14_init_state(params)

state = struct();

if params.IC_type == 1
   % Zero initial condition
   
   state.u_1 = 0;
   state.u_2 = 0;
   state.v   = 0;
   state.w_u = 0;
    
end
            
end

