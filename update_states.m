function [initial_states,U0,Vg_dot,Vg_total] = update_states(sampling_time,initial_states,u,f,generate_disturbance)

     dis_inflow = generate_disturbance;
     st = initial_states;
     con = u(1,:)'; % Apply first control move to th model to update the states

     [k1,Vg_dot,Vg_total] = f(st,con,dis_inflow);
     k2 = f(st+k1.*sampling_time/2,con,dis_inflow);
     k3 = f(st+k2.*sampling_time/2,con,dis_inflow);
     k4 = f(st+k3.*sampling_time,con,dis_inflow);
     st = st +sampling_time/6*(k1+2.*k2+2.*k3+k4);

     initial_states = full(st);
     
     %->....Removing first control and repeating last control
     U0 = [u(2:size(u,1),:);u(size(u,1),:)];


end