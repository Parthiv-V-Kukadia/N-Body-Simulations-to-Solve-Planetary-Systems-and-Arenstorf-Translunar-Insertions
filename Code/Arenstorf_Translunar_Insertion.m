    %% Initialize Values

clc
clear all
%Earth %Moon %Other Satellites/ Capsules

mu = 0.012277471;
mu_t = 1 - mu;
x_i = 2.005-0.994;
y_i = 0;
u_i = 0;
v_i = -26.91/40*2.001585106;

Time = 2.895/4*17.0752166; %Total time elapsed
steps = 60000; %Total number of points

%%
x = zeros(steps, 2);
v = zeros(steps, 2);

x(1,:) = [x_i y_i];
vel(1,:) = [u_i v_i];

t = 0;
dt = Time/steps; %time step (s)


for k = 1:steps
    x_i1 = x(k, 1); %Get the position in m
    y_i1 = x(k, 2);
    u_i1 = vel(k, 1); %Get the vel in m/s
    v_i1 = vel(k, 2); %Get the vel in m/s

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

    %update components
    x_ip1 = x_ip1 + dt/6*(k1x+2*k2x+2*k3x+k4x);
    y_ip1 = y_ip1 + dt/6*(k1y+2*k2y+2*k3y+k4y);
    u_ip1 = u_ip1 + dt/6*(k1u+2*k2u+2*k3u+k4u);
    v_ip1 = v_ip1 + dt/6*(k1v+2*k2v+2*k3v+k4v);

    x(k+1, 1) = x_ip1;
    x(k+1, 2) = y_ip1;
    vel(k+1, 1) = u_ip1;
    vel(k+1, 2) = v_ip1;
end

%%
figure(27)
subplot(1,5,[1:2])
%plot circle
theta = linspace(0, 2*pi, 100);
x_circ = cos(theta);
y_circ = sin(theta);
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',7)
hold on
plot(x_circ*384390.1/1000, y_circ*384390.1/1000, 'k--','LineWidth',3)
hold on
plot(x(:,1)*384390.1/1000,x(:,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [30 13])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [2 2 24 13])
set(gcf, 'PaperPosition', [2 2 24 13])
xlim([-500 500])
ylim([-500 500])
set(gca, 'FontSize',14)
x_val = [-500 -250 0 250 500];
y_val = [-500 -250 0 250 500];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
ylabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
text(250,0,'$1$','FontSize',16,'interpreter','latex');
text(-230,-170,'$2$','FontSize',16,'interpreter','latex');
text(70,250,'$3$','FontSize',16,'interpreter','latex');
text(70,-250,'$4$','FontSize',16,'interpreter','latex');
text(-230,170,'$5$','FontSize',16,'interpreter','latex');
%%
subplot(1,5,[3:5])
%plot circle
theta = linspace(0, 2*pi, 500);
x_circ = cos(theta);
y_circ = sin(theta);
plot(0,0,'.','Color',[0.2, 0.5470, 0.710],'MarkerSize',40)
hold on
plot(358.9,-132,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(x_circ*382390.1/1000, y_circ*382390.1/1000, 'k--','LineWidth',3)
hold on
plot(x(1:4278,1)*384390.1/1000,x(1:4278,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(x(steps-4340:steps,1)*384390.1/1000,x(steps-4340:steps,2)*384390.1/1000,'-','Color',[0.4, 0.6470, 0.410],'LineWidth',3)
hold on
plot(358.9,132,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(381.4,0,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(358.9,-132,'o','Color','k','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','w')
set(gcf, 'PaperUnits', 'centimeters')

set(gcf, 'Position', [2 2 36 13])
set(gcf, 'PaperPositionMode', 'auto')
xlim([-50 450])
ylim([-160 160])
set(gca, 'FontSize',14)
x_val = [0 100 200 300 400];
y_val = [-160 -80 0 80 160];
set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
xlabel('Distance ($\times1000$ km)','interpreter','latex','FontSize',16)
legend({'Earth','Moon Positions','Moon Orbit','Probe Orbit'},'Location','northwest','FontSize',14,'interpreter','latex')
annotation('arrow',[0.692 0.696], [0.774 0.777],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.696 0.692], [0.775-0.515 0.777-0.515],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.502-.15 0.506-.15], [0.60 0.60],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.2735 0.2735], [0.34 0.335],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.3436 0.3440], [0.34 0.335],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.2549 0.2536], [0.605 0.61],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.1761 0.1746], [0.575 0.58],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
annotation('arrow',[0.489 0.4896], [0.737 0.737],'HeadLength',15,'HeadWidth',15,'Color',[0.4, 0.6470, 0.410])
saveas(gcf,'TLI.png')
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