function [obj_total,g] = compute_both(U,H,P,Np,sampling_time,x_LRV_min,w_R,w_delU,w_U,f,n_states,n_u_prev)
    
    addpath('C:\aaa._MY_FILES\aUSN\4th semester\casadi-windows-matlabR2016a-v3.5.5');
    import casadi.*
   
    %->....Disturbance vector is assigned to P vector, these values stored 
    %in P vector are changed in each iteration fromt the main for loop
    
    %->....Objective function vector for each of the ensemble
    obj_array = SX.zeros(1,50);
    %->....Empty Constraints vector,later we will add the evaluated 
    %constraints in each step of prediction horizon through for loop
    g  = []; 
    
    %->....st variable initialization with the optimal states stored on H
    st = H(:,1); %Initial states
    
    %->....Initial condition constraints, the difference is appended in 
    %each iteration, the concept of multiple shooting adds this equality
    %constraint
    g = [g; st - P(1:n_states,1)];
    
   u_prev = P(3:4,1);
   ensemble_iteration = 0;
   
   for j = 1:50
        
        obj = 0;
        
        Vi_disturbance = P(n_states + n_u_prev + (j-1)*Np + 1:n_states + n_u_prev + j*Np,1);
   
        ensemble_iteration = ensemble_iteration + 1;
        
        sprintf('We are using %d numbered ensemble in the MPC',ensemble_iteration)
        
        U_extended = [u_prev,U];
        
        del_con = U_extended(:,2:end)-U_extended(:,1:end-1);
        
        for k = 1:Np
        
             st = H(:,k);
             con = U(:,k); 


            inflow = Vi_disturbance(k,1);

            %->....level of Merkebekk to be maximized in the objective function
            x_M1 = H(1,k) + x_LRV_min ; 
            obj = obj + (-w_R*x_M1^2 + del_con(:,k)'*...
                                       w_delU*del_con(:,k) + con'*w_U*con);


            st_next = H(:,k+1); 

            % Integrating the models with explicit runge kutta method
             k1 = f(st,con,inflow);
             k2 = f(st+k1.*sampling_time/2,con,inflow);
             k3 = f(st+k2.*sampling_time/2,con,inflow);
             k4 = f(st+k3.*sampling_time,con,inflow);

             st_next_predicted = st +sampling_time/6*(k1+2.*k2+2.*k3+k4);
             

             
              g = [g;
                     st_next - st_next_predicted % x(k+1) = f(x(k),u(k),d(k)) Equality constraints
                     ];               
        end
        obj_array(j) = obj; 
   end
   
    obj_total = sum(obj_array);
    
    end

   