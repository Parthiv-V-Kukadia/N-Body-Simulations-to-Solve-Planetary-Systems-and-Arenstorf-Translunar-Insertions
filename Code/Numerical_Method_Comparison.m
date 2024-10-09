    %% Initialize Values

clc
clear all
%Earth %Moon %Other Satellites/ Capsules

mu = 0.012277471;
mu_t = 1 - mu;
x_i = 0.994;
y_i = 0;
u_i = 0;
v_i = 1.00005*-2.001585106;

Time = 17.0752166; %Total time elapsed
steps = 15000; %Total number of points
%% Euler
x_Eu = zeros(steps, 2);
vel_Eu = zeros(steps, 2);

x_Eu(1,:) = [x_i y_i];
vel_Eu(1,:) = [u_i v_i];

t = 0;
dt = Time/steps; %time step (s)


for k = 1:steps
    x_i1 = x_Eu(k, 1); %Get the position in m
    y_i1 = x_Eu(k, 2);
    u_i1 = vel_Eu(k, 1); %Get the vel in m/s
    v_i1 = vel_Eu(k, 2); %Get the vel in m/s

    x_ip1 = x_i1;
    y_ip1 = y_i1;
    u_ip1 = u_i1;
    v_ip1 = v_i1;

    x_ip1 = x_ip1 + dt*u_ip1;
    y_ip1 = y_ip1 + dt*v_ip1;
    u_ip1 = u_ip1 + dt*du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);
    v_ip1 = v_ip1 + dt*dv_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);

    x_Eu(k+1, 1) = x_ip1;
    x_Eu(k+1, 2) = y_ip1;
    vel_Eu(k+1, 1) = u_ip1;
    vel_Eu(k+1, 2) = v_ip1;
end
%% Euler Plot
figure(7)
subplot(2,2,1)
%plot circle
theta = linspace(0, 2*pi, 100);
x_circ = cos(theta);
y_circ = sin(theta);
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_Eu(:,1)*384390.1/1000,x_Eu(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
xlim([-500 500])
ylim([-500 500])
legend({'Earth','Moon (i)','Probe (i)','Moon Orbit','Probe Orbit'},'Location','northeast','FontSize',11,'interpreter','latex')
title('Forward Euler','fontsize',16,'interpreter','latex')
annotation('arrow',[0.382 0.378], [0.741 0.741],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.249 0.246], [0.6964 0.691],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.339 0.34], [0.807 0.807],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
%% Heun
x_He = zeros(steps, 2);
vel_He = zeros(steps, 2);

x_He(1,:) = [x_i y_i];
vel_He(1,:) = [u_i v_i];

t = 0;
dt = Time/steps; %time step (s)


for k = 1:steps
    x_i1 = x_He(k, 1); %Get the position in m
    y_i1 = x_He(k, 2);
    u_i1 = vel_He(k, 1); %Get the vel in m/s
    v_i1 = vel_He(k, 2); %Get the vel in m/s

    xbar_i1 = x_i1 + dt*dx_dt(u_i1);
    ybar_i1 = y_i1 + dt*dy_dt(v_i1);
    ubar_i1 = u_i1 + dt*du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);
    vbar_i1 = v_i1 + dt*du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);
    
    x_ip1 = x_i1;
    y_ip1 = y_i1;
    u_ip1 = u_i1;
    v_ip1 = v_i1;

    x_ip1 = x_ip1 + (dt/2)*(dx_dt(u_i1) + dx_dt(ubar_i1));
    y_ip1 = y_ip1 + (dt/2)*(dy_dt(v_i1) + dy_dt(vbar_i1));
    u_ip1 = u_ip1 + (dt/2)*(du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t)...
        + du_dt(xbar_i1, ubar_i1, ybar_i1, vbar_i1, mu, mu_t));
    v_ip1 = v_ip1 + (dt/2)*(dv_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t)...
        + dv_dt(xbar_i1, ubar_i1, ybar_i1, vbar_i1, mu, mu_t));

    x_He(k+1, 1) = x_ip1;
    x_He(k+1, 2) = y_ip1;
    vel_He(k+1, 1) = u_ip1;
    vel_He(k+1, 2) = v_ip1;
end

%% Heun Plot
subplot(2,2,2)
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_He(1:11800,1)*384390.1/1000, x_He(1:11800,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [8 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
title("Heun's Method",'fontsize',16,'interpreter','latex')
annotation('arrow',[0.826 0.822], [0.7415 0.7415],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.812 0.808], [0.690 0.6930],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.756 0.7555], [0.820 0.8234],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.612-0.025 0.608-0.025], [0.710 0.7136],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.678-0.025 0.675-0.025], [0.610 0.6136],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.782 0.786], [0.7185 0.7185],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
%% KDK

f = 1.18;
Time = 17.0752166*f;

x_KDK = zeros(steps*f, 2);
vel_KDK = zeros(steps*f, 2);

x_KDK(1,:) = [x_i y_i];
vel_KDK(1,:) = [u_i v_i];

t = 0;
dt = Time/steps/f; %time step (s)


for k = 1:steps*f
    x_i1 = x_KDK(k, 1); %Get the position in m
    y_i1 = x_KDK(k, 2);
    u_i1 = vel_KDK(k, 1); %Get the vel in m/s
    v_i1 = vel_KDK(k, 2); %Get the vel in m/s

    u_i12 = u_i1 + du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t)*dt/2;
    v_i12 = v_i1 + dv_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t)*dt/2;
    
    x_ip1 = x_i1;
    y_ip1 = y_i1;
    u_ip1 = u_i12;
    v_ip1 = v_i12;

    x_ip1 = x_ip1 + dt*u_i12;
    y_ip1 = y_ip1 + dt*v_i12;
    u_ip1 = u_ip1 + du_dt(x_ip1, u_ip1, y_ip1, v_ip1, mu, mu_t)*dt/2;
    v_ip1 = v_ip1 + dv_dt(x_ip1, u_ip1, y_ip1, v_ip1, mu, mu_t)*dt/2;

    x_KDK(k+1, 1) = x_ip1;
    x_KDK(k+1, 2) = y_ip1;
    vel_KDK(k+1, 1) = u_ip1;
    vel_KDK(k+1, 2) = v_ip1;
end

%% KDK Plot
subplot(2,2,3)
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_KDK(:,1)*384390.1/1000, x_KDK(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
title('KDK','fontsize',16,'interpreter','latex')
annotation('arrow',[0.3105 0.305], [0.79-0.455 0.795-0.457],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.2403 0.237], [0.308 0.303],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.2505 0.26], [0.243 0.234],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.335 0.34], [0.242 0.251],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow', [0.258 0.2495], [0.183 0.173],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
%% RK4

Time = 17.0752166;

x_RK4 = zeros(steps, 2);
vel_RK4 = zeros(steps, 2);

x_RK4(1,:) = [x_i y_i];
vel_RK4(1,:) = [u_i v_i];

t = 0;
dt = Time/steps; %time step (s)


for k = 1:steps
    x_i1 = x_RK4(k, 1); %Get the position in m
    y_i1 = x_RK4(k, 2);
    u_i1 = vel_RK4(k, 1); %Get the vel in m/s
    v_i1 = vel_RK4(k, 2); %Get the vel in m/s

    x_ip1 = x_i1;
    y_ip1 = y_i1;
    u_ip1 = u_i1;
    v_ip1 = v_i1;

    % first round
    k1x = dx_dt(u_i1);
    k1y = dy_dt(v_i1);
    k1u = du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);
    k1v = dv_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);

    % second round
    k2x = dx_dt(u_i1+dt*k1u/2);
    k2y = dy_dt(v_i1+dt*k1v/2);
    k2v = dv_dt(x_i1+dt*k1x/2, u_i1+dt*k1u/2, y_i1+dt*k1y/2, v_i1+dt*k1v/2, mu, mu_t);
    k2u = du_dt(x_i1+dt*k1x/2, u_i1+dt*k1u/2, y_i1+dt*k1y/2, v_i1+dt*k1v/2, mu, mu_t);

    % third round
    k3x = dx_dt(u_i1+dt*k2u/2);
    k3y = dy_dt(v_i1+dt*k2v/2);
    k3v = dv_dt(x_i1+dt*k2x/2, u_i1+dt*k2u/2, y_i1+dt*k2y/2, v_i1+dt*k2v/2, mu, mu_t);
    k3u = du_dt(x_i1+dt*k2x/2, u_i1+dt*k2u/2, y_i1+dt*k2y/2, v_i1+dt*k2v/2, mu, mu_t);

    % fourth round
    k4x = dx_dt(u_i1+dt*k3u);
    k4y = dy_dt(v_i1+dt*k3v);
    k4v = dv_dt(x_i1+dt*k3x, u_i1+dt*k3u, y_i1+dt*k3y, v_i1+dt*k3v, mu, mu_t);
    k4u = du_dt(x_i1+dt*k3x, u_i1+dt*k3u, y_i1+dt*k3y, v_i1+dt*k3v, mu, mu_t);

    x_ip1 = x_ip1 + dt/6*(k1x+2*k2x+2*k3x+k4x);
    y_ip1 = y_ip1 + dt/6*(k1y+2*k2y+2*k3y+k4y);
    u_ip1 = u_ip1 + dt/6*(k1u+2*k2u+2*k3u+k4u);
    v_ip1 = v_ip1 + dt/6*(k1v+2*k2v+2*k3v+k4v);

    x_RK4(k+1, 1) = x_ip1;
    x_RK4(k+1, 2) = y_ip1;
    vel_RK4(k+1, 1) = u_ip1;
    vel_RK4(k+1, 2) = v_ip1;
end
%% RK4 Plot
subplot(2,2,4)
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_RK4(:,1)*384390.1/1000,x_RK4(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
title('RK4','fontsize',16,'interpreter','latex')
annotation('arrow',[0.7905 0.785], [0.79-0.474 0.795-0.474],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.795 0.8005], [0.724-0.474 0.729-0.474],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.684 0.680], [0.796-0.474 0.790-0.474],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.690 0.694], [0.703-0.474 0.697-0.474],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,'Arenstorf.png')

%% Stability Plot- RK4
figure(987)
subplot(2,2,1)
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_RK4(:,1)*384390.1/1000,x_RK4(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
title('1 Period','fontsize',16,'interpreter','latex')

subplot(2,2,2)
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_RK4(:,1)*384390.1/1000,x_RK4(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
plot(x_RK4(:,1)*384390.1/1000*1.02,x_RK4(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
title('2 Periods','fontsize',16,'interpreter','latex')


%% RK4 Iteration

Time = 2*17.0752166;

x_RK4 = zeros(steps, 2);
vel_RK4 = zeros(steps, 2);

x_RK4(1,:) = [x_i y_i];
vel_RK4(1,:) = [u_i v_i];

t = 0;
dt = Time/2/steps; %time step (s)


for k = 1:2*steps
    x_i1 = x_RK4(k, 1); %Get the position in m
    y_i1 = x_RK4(k, 2);
    u_i1 = vel_RK4(k, 1); %Get the vel in m/s
    v_i1 = vel_RK4(k, 2); %Get the vel in m/s

    x_ip1 = x_i1;
    y_ip1 = y_i1;
    u_ip1 = u_i1;
    v_ip1 = v_i1;

    % first round
    k1x = dx_dt(u_i1);
    k1y = dy_dt(v_i1);
    k1u = du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);
    k1v = dv_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);

    % second round
    k2x = dx_dt(u_i1+dt*k1u/2);
    k2y = dy_dt(v_i1+dt*k1v/2);
    k2v = dv_dt(x_i1+dt*k1x/2, u_i1+dt*k1u/2, y_i1+dt*k1y/2, v_i1+dt*k1v/2, mu, mu_t);
    k2u = du_dt(x_i1+dt*k1x/2, u_i1+dt*k1u/2, y_i1+dt*k1y/2, v_i1+dt*k1v/2, mu, mu_t);

    % third round
    k3x = dx_dt(u_i1+dt*k2u/2);
    k3y = dy_dt(v_i1+dt*k2v/2);
    k3v = dv_dt(x_i1+dt*k2x/2, u_i1+dt*k2u/2, y_i1+dt*k2y/2, v_i1+dt*k2v/2, mu, mu_t);
    k3u = du_dt(x_i1+dt*k2x/2, u_i1+dt*k2u/2, y_i1+dt*k2y/2, v_i1+dt*k2v/2, mu, mu_t);

    % fourth round
    k4x = dx_dt(u_i1+dt*k3u);
    k4y = dy_dt(v_i1+dt*k3v);
    k4v = dv_dt(x_i1+dt*k3x, u_i1+dt*k3u, y_i1+dt*k3y, v_i1+dt*k3v, mu, mu_t);
    k4u = du_dt(x_i1+dt*k3x, u_i1+dt*k3u, y_i1+dt*k3y, v_i1+dt*k3v, mu, mu_t);

    x_ip1 = x_ip1 + dt/6*(k1x+2*k2x+2*k3x+k4x);
    y_ip1 = y_ip1 + dt/6*(k1y+2*k2y+2*k3y+k4y);
    u_ip1 = u_ip1 + dt/6*(k1u+2*k2u+2*k3u+k4u);
    v_ip1 = v_ip1 + dt/6*(k1v+2*k2v+2*k3v+k4v);

    x_RK4(k+1, 1) = x_ip1;
    x_RK4(k+1, 2) = y_ip1;
    vel_RK4(k+1, 1) = u_ip1;
    vel_RK4(k+1, 2) = v_ip1;
end
%%

subplot(2,2,3)
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_RK4(:,1)*384390.1/1000,x_RK4(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
title('3 Periods','fontsize',16,'interpreter','latex')


%% RK4 Iteration Mega
f =3;
Time = f*17.0752166;

x_RK4 = zeros(steps, 2);
vel_RK4 = zeros(steps, 2);

x_RK4(1,:) = [x_i y_i];
vel_RK4(1,:) = [u_i v_i];

t = 0;
dt = Time/f/steps; %time step (s)


for k = 1:f*steps
    x_i1 = x_RK4(k, 1); %Get the position in m
    y_i1 = x_RK4(k, 2);
    u_i1 = vel_RK4(k, 1); %Get the vel in m/s
    v_i1 = vel_RK4(k, 2); %Get the vel in m/s

    x_ip1 = x_i1;
    y_ip1 = y_i1;
    u_ip1 = u_i1;
    v_ip1 = v_i1;

    % first round
    k1x = dx_dt(u_i1);
    k1y = dy_dt(v_i1);
    k1u = du_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);
    k1v = dv_dt(x_i1, u_i1, y_i1, v_i1, mu, mu_t);

    % second round
    k2x = dx_dt(u_i1+dt*k1u/2);
    k2y = dy_dt(v_i1+dt*k1v/2);
    k2v = dv_dt(x_i1+dt*k1x/2, u_i1+dt*k1u/2, y_i1+dt*k1y/2, v_i1+dt*k1v/2, mu, mu_t);
    k2u = du_dt(x_i1+dt*k1x/2, u_i1+dt*k1u/2, y_i1+dt*k1y/2, v_i1+dt*k1v/2, mu, mu_t);

    % third round
    k3x = dx_dt(u_i1+dt*k2u/2);
    k3y = dy_dt(v_i1+dt*k2v/2);
    k3v = dv_dt(x_i1+dt*k2x/2, u_i1+dt*k2u/2, y_i1+dt*k2y/2, v_i1+dt*k2v/2, mu, mu_t);
    k3u = du_dt(x_i1+dt*k2x/2, u_i1+dt*k2u/2, y_i1+dt*k2y/2, v_i1+dt*k2v/2, mu, mu_t);

    % fourth round
    k4x = dx_dt(u_i1+dt*k3u);
    k4y = dy_dt(v_i1+dt*k3v);
    k4v = dv_dt(x_i1+dt*k3x, u_i1+dt*k3u, y_i1+dt*k3y, v_i1+dt*k3v, mu, mu_t);
    k4u = du_dt(x_i1+dt*k3x, u_i1+dt*k3u, y_i1+dt*k3y, v_i1+dt*k3v, mu, mu_t);

    x_ip1 = x_ip1 + dt/6*(k1x+2*k2x+2*k3x+k4x);
    y_ip1 = y_ip1 + dt/6*(k1y+2*k2y+2*k3y+k4y);
    u_ip1 = u_ip1 + dt/6*(k1u+2*k2u+2*k3u+k4u);
    v_ip1 = v_ip1 + dt/6*(k1v+2*k2v+2*k3v+k4v);

    x_RK4(k+1, 1) = x_ip1;
    x_RK4(k+1, 2) = y_ip1;
    vel_RK4(k+1, 1) = u_ip1;
    vel_RK4(k+1, 2) = v_ip1;
end
%%
subplot(2,2,4)
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x_RK4(:,1)*384390.1/1000,x_RK4(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(384390.1/1000,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(0.994*384390.1/1000,0,'r.','MarkerSize',30)
set(gca, 'FontSize',14)
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500 750];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
title('4 Periods','fontsize',16,'interpreter','latex')

set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,'Aren_Stability.png')
%%
function ans = du_dt(x, x_t, y, y_, mu, mu_t)
    D1 = ((x+mu)^2+y^2)^(3/2);
    D2 = ((x-mu_t)^2+y^2)^(3/2);
    ans = x + 2*y_ - mu_t*(x+mu)/D1 - mu*(x-mu_t)/D2;
end
%%
function ans = dv_dt(x, x_, y, y_, mu, mu_t)
    D1 = ((x+mu)^2+y^2)^(3/2);
    D2 = ((x-mu_t)^2+y^2)^(3/2);
    ans = y - 2*x_ - mu_t*y/D1 - mu*y/D2;
end
%%
function ans = dx_dt(u)
    ans = u;
end
%%
function ans = dy_dt(v)
    ans = v;
end