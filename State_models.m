function [h_dot, Vg_dot,Vg_dot_total]  = State_models(h,hg,Vi)

    %->....Parameters declaration
    alpha = 0.05; % Fraction of surface area in compartment 2
    beta = 0.02; % Fraction of inflow to compartment 2
    
    C_d = 0.7; % Discharge coefficient, Dalsfos gate
    w = [11.6,11]; % Width of the flood gates 1 and 2
    g = 9.8; % Acceleration due to gravity
    K12 = 800; % Inter compartmental flow coefficient
    
  
    
    %->....Disturbance assignment to a variable
    Vi_dot = Vi;
    %->....Volumetric flow through the turbine [m3/sec]
    Vt_dot = 40; 
    
    
    %->...Intermediate equations: Algebric equations
    sectional_area = @(h) max((28e6*1.1*h^(1/10)),1e3);
    V12_dot = K12*(h(1)-h(2))*sqrt(abs(h(1)-h(2))); % Intercompartmental flow 
    %->....Vg_dot = C_d.*w.*min(hg,h(2))'*sqrt(2*g*max(h(2),0));
    %->....flow through the flood gates
    Vg_dot = [C_d*w(1)*min(hg(1),h(2))*sqrt(2*g*max(h(2),0));
              C_d*w(2)*min(hg(2),h(2))*sqrt(2*g*max(h(2),0))
              ];
    %->....Total flow through both gates
    Vg_dot_total = sum(Vg_dot);
    
    % State equations
    
    h_dot = [((1-beta)*Vi_dot-V12_dot)/((1-alpha)*sectional_area(h(1)));
        
            (beta*Vi_dot-Vt_dot-Vg_dot_total+V12_dot)/(alpha*sectional_area(h(2))) ];

        
end



