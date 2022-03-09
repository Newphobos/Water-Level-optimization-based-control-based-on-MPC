function [obj,g] = compute_both(U,H,P,Np,sampling_time,x_LRV_min,w_R,w_delU,w_U,f)
    

    %->....Objective function initialization
    obj = 0;
    %->....Empty Constraints vector,later we will add the evaluated 
    %constraints in each step of prediction horizon through for loop
    g  = []; 
    
    %->....st variable initialization with the optimal states stored on H
    st = H(:,1); %Initial states
    
    %->....Initial condition constraints, the difference is appended in 
    %each iteration, the concept of multiple shooting adds this equality
    %constraint
    g = [g; st - P(1:2,1)]; 
    
    %->....Disturbance vector is assigned to P vector, these values stored 
    %in P vector are changed in each iteration fromt the main for loop
    Vi_disturbance = P(3:end,1);
    
    
 for k = 1:Np
        
        st = H(:,k);
        con = U(:,k); 
        con_kplus1 = U(:,k+1);                             
        del_con = con_kplus1 - con;
        
        inflow = Vi_disturbance(k,1);
        
        %->....level of Merkebekk to be maximized in the objective function
        x_M1 = H(1,k) + x_LRV_min ; 
        obj = obj + (-w_R*x_M1^2 + (con_kplus1-con)'*...
                                   w_delU*(con_kplus1-con) + con'*w_U*con);
 

        st_next = H(:,k+1); 
        
        % Integrating the models with explicit runge kutta method
         k1 = f(st,con,inflow);
         k2 = f(st+k1.*sampling_time/2,con,inflow);
         k3 = f(st+k2.*sampling_time/2,con,inflow);
         k4 = f(st+k3.*sampling_time,con,inflow);
         
         st_next_predicted = st +sampling_time/6*(k1+2.*k2+2.*k3+k4);
          
         g = [g;
             st_next - st_next_predicted; % x(k+1) = f(x(k),u(k)) Equality constraints
             del_con % Ineqality constraint for change in control signal
             ];              
 end

end