function [] = plotting_results(x,u_eachstep,Vgate_each,Vgate_total,...
                                generate_disturbance,xLRV,...
                            xHRV,timesteps,time_eachloop,...
                            Average_looprun_time,rl,ru)

%     close all
    %->....Plot results showing the level changes x_M and x_D
    figure(1)
    subplot(4,1,1,'align');
    plot5 = plot(1:1:timesteps+1,x,'LineWidth',1.05); 
    hold on;
    plot6 = plot(1:1:timesteps,xLRV(1:timesteps,1),'-r',...
                 1:1:timesteps,xHRV(1:timesteps,1),'-r','LineWidth',1.05);
    hold on;
    plotn = plot(1:1:timesteps,rl(1:timesteps,1),'--k',...
                 1:1:timesteps,ru(1:timesteps,1),'--k','LineWidth',1.05);
    combined_plot = [plot5;plot6;plotn];
    legend(combined_plot,'x_M','x_D','xLRV','xHRV','r_l','r_u');
    xlabel('Time(hours)');
    ylabel('h[m]');
    title('Initial Levels = [58.25;58.05]')
    grid on;

    
    %->....Plot of the control signal,hg in each time step
    subplot(4,1,2,'align');
    plot(1:1:timesteps,u_eachstep,'LineWidth',1);
%     ylim([0,600]);
    legend('hg1','hg2');
    xlabel('Time(hours)');
    ylabel('hg[m]');
%     title('Gate opening')
    grid on;
    
    %->....Plot of the volumetric flow through the gates
    subplot(4,1,3,'align');
    plota = plot(1:1:timesteps,Vgate_each,'LineWidth',0.95); 
    hold on;
    plotb = plot(1:1:timesteps,Vgate_total,'-r','LineWidth',0.95);
    hold on;
    combined_plot = [plota;plotb];
    legend(combined_plot,'Vg1','Vg2','Total Vg');
    xlabel('Time(hours)');
    ylabel('flow[m3/sec]');
%     ylim([0,200]);
    grid on;
%     title('Flow through gates')
    
    %->....Plot disturbance flow, Vi
    subplot(4,1,4,'align')
    plot(1:1:timesteps-1,generate_disturbance(1:1:timesteps-1),'LineWidth',0.95)
    xlabel('Time(hours)');
    ylabel('Vi[m3/sec]');
    legend('Vi');
    grid on;
%     title('inflow into the lake')

    %->....Time elapsed in each loop execution
    figure(2)
    plot9 = plot(time_eachloop);
    hold on;
    plot8 = plot(ones(timesteps,1)*Average_looprun_time,'-g');
    plots = [plot9;plot8];
    legend(plots,'time each loop','Average time each loop');
    xlabel('Number of loop');
    ylabel('time[sec]')
    title('Time elapsed in each iteration');
    
    figure(3)
    plot5 = plot(1:1:timesteps+1,x,'LineWidth',1.05); 
    hold on;
    plot6 = plot(1:1:timesteps,xLRV(1:timesteps,1),'--r',...
                 1:1:timesteps,xHRV(1:timesteps,1),'--r','LineWidth',1.05);
    hold on;
   
    combined_plot = [plot5;plot6];
    legend(combined_plot,'x_M','x_D','xLRV','xHRV');
    xlabel('Time(hours)');
    ylabel('h[m]');
    title('Initial Levels = [58.25;58.05]')
    grid on;
    

end