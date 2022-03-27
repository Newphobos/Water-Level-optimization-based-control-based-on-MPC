

clc
clear
% close all

addpath('C:\aaa._MY_FILES\aUSN\4th semester\casadi-windows-matlabR2016a-v3.5.5');
import casadi.*


%%....................Parameters settings::Start.........................%%

%->....Gate opening requirements
hg_min = 0; % [m], minimum gate opening
hg_max = 5.6;%[m], maximum gate opening

%->....Time setting for simulations
sampling_time = 3600;% [s],Sampling time of 1 hour
sim_time = 30*24*3600;% [s],Simulation time of 30 days
timespan = 0:sampling_time:sim_time;
timesteps = length(timespan);

%->....Length of the prediction horizon
Np = 13*24; % 13 days of prediction horizon

%->....Level constraints from April 15 to May 15
x_LRV1 = 55.75; %[m] from april 15 to april 30    
x_LRV2 = 58.85; % [m] from may 1st to may 15
x_HRV1 = 60.35; %[m] from april 15 to april 30
x_HRV2 = 59.85; % [m] from may 1st to may 15

%->....Creating array for constraints limit and also needed for plotting
xLRV = zeros(timesteps,1);
xLRV(1:timesteps/2) = x_LRV1;
xLRV(timesteps/2:end) = x_LRV2;

xHRV = zeros(timesteps,1);
xHRV(1:timesteps/2) = x_HRV1;
xHRV(timesteps/2:end) = x_HRV2;

%->....Add extra number of values for simulating until the end of...
%->simulation time
xLRV = [xLRV;ones(Np+1,1).*x_LRV2];
xHRV = [xHRV;ones(Np+1,1).*x_HRV2];

%->....Finding minimum value of XLRV, added to states to find x_M and x_D
x_LRV_min = min(xLRV); %[m] 


%->.... Nominal values for MPC cost function parameters
X_R = 0.75; % Location of rl_i in interval [xLRV,xHRV]
delta_HRV = 0.05; %[m] Safety margin for upper reference region boundary

%Selection of reference region boundries
 rl = (1-X_R).*xLRV + X_R.* xHRV;
 ru = xHRV - delta_HRV;

%->.... Nominal values for MPC cost function parameters
X_R = 0.75; % Location of rl_i in interval [xLRV,xHRV]
delta_HRV = 0.05; %[m] Safety margin for upper reference region boundary

%Weight matrices
w_R = 10; %Weight on cost outside reference region
w_delU = 1; %Weight on change of flood gate opening
w_U = 0.1 ; % Weight on flood gate opening
    
% Vectors for storing volumetric flow through gates
Vgate_each = zeros(timesteps,2);
Vgate_total = zeros(timesteps,1);

%%...........Parameters settings::End................................%%

%%..........Symbolic representation for CasADi::Start................%%

%->....States representation,two states h1 and h2
h = SX.sym('h',2,1);
st = h ;
n_states = length(st);

%->....Control action representation,hg1 and hg2
hg = SX.sym('hg',2,1);
controls = hg;
n_controls = length(controls);
n_u_prev = n_controls;

%->....Process disturbance, volumetric inflow
Vi_dot = SX.sym('Vi_dot');
n_disturbance = length(Vi_dot);

%->....Process model, rhs term and the output from the models 
[hdot,Vg_dot,Vg_dot_total] = State_models(st,controls,Vi_dot); 

%->.... Non Linear mapping function; inputs to outputs
f = Function('f',{st,controls,Vi_dot},{hdot,Vg_dot,Vg_dot_total}); 

%->....Decision variables(control input) and states vector for the whole
%->prediction horizon
U = SX.sym('U',n_controls,Np); % Control inputs matrix
H = SX.sym('H',n_states,(Np + 1)); % States matrix
 


%->.... Parameter vector to store initial states and disturbance
P = SX.sym('P',n_states + n_u_prev + 50*Np,1); % 2 initial states and Np number of 
%->disturbance for the whole prediction horiozon

%%..........Symbolic representation for CasADi::END................%%

%%..........Evalution of objective and constraint::Start................%%

%->....Compute objective and constraints over the prediction horizon
[obj,g] = compute_both(U,H,P,Np,sampling_time,x_LRV_min,w_R,w_delU,w_U,f,n_states,n_u_prev);
%%..........Evalution of objective and constraint::END................%%

%%..........IPOPT optimizer setup::Start................%%

%->....Define decision/optimization variables
OPT_variables = [reshape(H,2*(Np+1),1);reshape(U,2*Np,1)];
 
%->....Define NLP problem object
nlp_prob = struct('f',obj,'x',OPT_variables,'g',g,'p',P);

%->....Options set up for optimizer
opts = struct;
opts.ipopt.max_iter = 1000;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-6;
opts.ipopt.acceptable_obj_change_tol = 1e-6;
opts.ipopt.print_level = 0; % 3,5(default) upto 13
opts.ipopt.fixed_variable_treatment = 'make_constraint';
opts.ipopt.warm_start_init_point = 'yes';

solver = nlpsol('solver','ipopt',nlp_prob,opts); 
  
%%..........IPOPT optimizer setup::END................%%

%%..........Bounds setup::Start................%%
args = struct;
%->....equality constrainst matrix g(equality constraint) and delta_u 
% -0.01 <= delta_u <=0.01
% args.lbg = repelem([0;-0.01;0],[2*(Np+1),2*Np,49*(2*Np)]);
% args.ubg = repelem([0;0.01;0],[2*(Np+1),2*Np,49*(2*Np)]);

args.lbg = repelem(0,2*(50*Np+1));
args.ubg = repelem(0,2*(50*Np+1));


%->....Bounds on gate opening, U
args.lbx(2*(Np+1)+1:1:2*(Np+1)+ 2*Np,1) = 0; % Lower limit of gate opening
args.ubx(2*(Np+1)+1:1:2*(Np+1)+ 2*Np,1) = hg_max;
%Upper bound of gate opening
%..........Bound setup::Partial END................%%
%.........Partial END because the bounds on states are changed in each
%%time step, therefore the bounds are supplied fromt the main loop.


%..........The simulation loop starts from here........................%%
% % --------------------------------------------------------------------------
% % --------------------------------------------------------------------------

%->.... Initialization of control input over the whole prediction horizon
initial_control = [0;0];
U0 = repmat(initial_control,1,Np)';
u_prev = [0;0];

%->....Initialization of states over the whole prediction horizon
%->....Give initial states between 0 and 4.6, this limit is because it will
%help to bound x_M and x_D inbetween xLRV and xHRV
initial_states = [2.0;1.95]; 
H0 =  repmat(initial_states,1,Np+1)'; 


%->....state_history matrix contains the evolution of states over the whole
%simulation time, these states values are obtained from the first optimized
%control input in each time step
state_history(:,1) = initial_states; 

%->....u_eachstep contains the changes in control signal over the whole 
%simulation time period, these values are first optimal move stored in each
%time step
u_eachstep = zeros(2,timesteps);

%->....array for storing each loop execution time

time_eachloop = zeros(timesteps,1);

%->.... Create a large vector storing the value for the disturbance over
%the simulation period
disturbance = select_disturbance();
%-> ....Repeat each element of inflow disturbance(Here,24 times), because
%prection of inflow from the company are provided after 24 hours, i.e next
%day
generate_disturbance = repelem(disturbance,24).*3;

%->.... Start MPC with with initialization of iterations count
no_iterations = 0;
%->....For optimal solution trajectory starting each timestep
% optimal_traj = []; % Contains the optimal solution trajectory



all_ensemble = access_ensemble();

days_counter = 1;
%->.....Start the main loop
for i=1:timesteps
     
    tic;
    %->....Parametrization vector,P and initial vector args.x0 are changed
    %in each timestep
%     args.p = [initial_states;generate_disturbance(i:i+Np-1,1)];
    
%     args.p = [initial_states;u_prev;all_ensemble(:,days_counter)];
    
    if  mod(i,24) == 0
        days_counter = days_counter + 1 ;
        args.p = [initial_states;u_prev;all_ensemble(:,days_counter)];
        
    else 
        
        args.p = [initial_states;u_prev;all_ensemble(:,days_counter)];
        
    end
    
    args.x0 = [reshape(H0',2*(Np+1),1);reshape(U0',2*Np,1)];

    %->....Bounds on h1 
    % x_M E [xLRV,xHRV]
    % xLRV <= x_M <= xHRV
    %xLRV <= h1+x_LRV_min <= xHRV
    %xLRV-x_LRV_min <= h1 <= xHRV-x_LRV_min
    
    args.lbx(1:2:2*(Np+1),1) = xLRV(i:i+Np,1) - x_LRV_min; 
    args.ubx(1:2:2*(Np+1),1) = xHRV(i:i+Np,1) - x_LRV_min;

%     %->....Bounds on h2
    args.lbx(2:2:2*(Np+1),1) = xLRV(i:i+Np,1)-x_LRV_min;
    args.ubx(2:2:2*(Np+1),1) = xHRV(i:i+Np,1)-x_LRV_min;

    %->.... Call IPOPT solver object in each iterations
    sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',...
                    args.lbg,'ubg',args.ubg,'p',args.p);
    
    %->....Extract the control signal which is after states variables in 
    % the optimization variable
    u = reshape(full(sol.x(2*(Np+1)+1:end))',2,Np)';
     %Get OPTIMAL solution trajectory  
     %optimal_traj(:,1:2,no_iterations+1) = 
     %reshape(full(sol.x(1:2*(Np+1)))',2,Np+1)';     
   
    %->....Append first move in eachstep
    u_eachstep(:,i) = u(1,:)';
    
    %->....Update the states by appling first control move
    [initial_states,U0,Vg_dot,Vg_total] = update_states(sampling_time,...
                                    initial_states,u,f,generate_disturbance(i));
    
    u_prev = u(1,:)';                         
                                
    %->....Storing gate flow in each time step
    Vgate_each(no_iterations+1,:) = full(Vg_dot);
    Vgate_total(no_iterations+1,:) = full(Vg_total);  

    % ->....Storing updated states in each time step
    state_history(:,no_iterations+2) = initial_states;
     
    %->....Get solution trajectory
    H0 = reshape(full(sol.x(1:2*(Np+1)))',2,Np+1)';
    
    %->....Shift trajectory to initialize the next step
    H0 = [H0(2:end,:);H0(end,:)];
    
    %->....Increase the number of iterations
    no_iterations = no_iterations + 1
    
    %->....Record time of completion of each loop 
    time_eachloop(i) = toc; 
end
 
%->....Finding average loop completion time
Average_looprun_time = sum(time_eachloop)/timesteps
Total_iteration = no_iterations;
total_time = Total_iteration*Average_looprun_time

%->.... x holds the information both x_M and x_D, the levels we are 
% interested on.
x = state_history + x_LRV_min;


% ->....Visualization function call for plotting all the results
plotting_results(x,u_eachstep,Vgate_each,Vgate_total,generate_disturbance,...
                  xLRV,xHRV,timesteps,time_eachloop,Average_looprun_time,rl,ru);



 