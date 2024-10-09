    %% Initialize Values
clc
clear all
%Earth %Moon %Other Satellites/ Capsules

N = 2; %number of bodies (<= 5)

Time = 900; %Total time elapsed (hours)
step = 1000*[1 10 100 1000]; %Total number of points

%Change up to the Nth Initial Body Position
x1_i = [-1.23e7 3.84e7]; %Initial 1st Body Position (m)
x2_i = [1.23e7 -3.84e7]; %Initial 2nd Bodu Position (m)

%Change up to the Nth Initial Body Velocity
v1_i = [-380 120]; %Initial 1st Body Velocity (m/s)
v2_i = [380 -120]; %Initial 2nd Body Velocity (m/s)

%Change up to the Nth Initial Body Mass
m_1 = 5.972e23; %1st Body Mass (kg)
m_2 = 5.972e23; %2nd Body Mass (kg)

%%
for alpha = 1:length(step)
    steps = step(alpha);
    x = zeros(steps, 4);
    v = zeros(steps, 4);
    
    t = 0;
    dt = Time*3600/steps; %time step (s)
    mass = [m_1, m_2];
    x(1,:) = [x1_i, x2_i];
    v(1,:) = [v1_i, v2_i];
    for k = 1:steps
        t = t + dt;
        for i = 1:N
            x_i = x(k, ((i-1)*2+1):i*2); %Get the position in m
            v_i = v(k, ((i-1)*2+1):i*2); %Get the vel in m/s

            x_ip1 = x_i;
            v_ip1 = v_i;
            const = 0;
            for j = 1:N
                if j ~= i
                    x_j = x(k,((j-1)*2+1):j*2);
                    
                    k1x = dX_dT(mass(1,j), x_i, x_j, v_i);
                    k1v = dV_dT(mass(1,j), x_i, x_j);
                    
                    k2x = dX_dT(mass(1,j), x_i+dt*k1x/2, x_j, v_i+dt*k1v/2);
                    k2v = dV_dT(mass(1,j), x_i + dt*k1x/2, x_j);
                    
                    k3x = dX_dT(mass(1,j), x_i+dt*k2x/2, x_j, v_i+dt*k2v/2);
                    k3v = dV_dT(mass(1,j), x_i + dt*k2x/2, x_j);
                    
                    k4x = dX_dT(mass(1,j), x_i+dt*k3x, x_j, v_i+dt*k3v);
                    k4v = dV_dT(mass(1,j), x_i + dt*k3x, x_j);

                    x_ip1 = x_ip1 + dt/6*(k1x+2*k2x+2*k3x+k4x);
                    v_ip1 = v_ip1 + dt/6*(k1v+2*k2v+2*k3v+k4v);
                end
            end
            x(k+1, ((i-1)*2+1):i*2) = x_ip1;
            v(k+1, ((i-1)*2+1):i*2) = v_ip1;
        end
    end
    figure(1568)
    subplot(2,2,alpha)
    plot(x(:,1)/10e5,x(:,2)/10e5,'LineWidth',3)
    hold on
    plot(x(:,3)/10e5,x(:,4)/10e5,'LineWidth',3)
    hold on
    xlabel('Distance ($\times100$ km)','interpreter','latex','fontsize',16)
    ylabel('Distance ($\times100$ km)','interpreter','latex','fontsize',16)
    if alpha==1
        set(gcf, 'PaperUnits', 'centimeters')
        set(gcf, 'PaperSize', [22 22])
        set(gcf, 'Units', 'centimeters' )
        set(gcf, 'Position', [20 0 22 22])
        set(gcf, 'PaperPosition', [0 0 22 22])
        xlim([-400 400])
        ylim([-400 400])
        x_val = [-400 -200 0 200 400];
        y_val = [-400 -200 0 200 400];
        set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
        set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
        title('$\Delta t = 3240$ s','interpreter','latex','fontsize',16)
    elseif alpha ==2
        x_val = [-100 -50 0 50 100];
        y_val = [-100 -50 0 50 100];
        set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
        set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
        title('$\Delta t = 324$ s','interpreter','latex','fontsize',16)
    elseif alpha ==3
        x_val = [-50 -25 0 25 50];
        y_val = [-50 -25 0 25 50];
        set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
        set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
        title('$\Delta t = 32.4$ s','interpreter','latex','fontsize',16)
    elseif alpha ==4
        x_val = [-50 -25 0 25 50];
        y_val = [-50 -25 0 25 50];
        set(gca,'xtick', x_val, 'xticklabel', num2str(x_val.'))
        set(gca,'ytick', y_val, 'yticklabel', num2str(y_val.'))
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16)
        title('$\Delta t = 3.24$ s','interpreter','latex','fontsize',16)
    end
end
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,'N=2_Simulation.png')
%%
function ans = dV_dT(m_j, r_i, r_j)
    G = 6.67430e-11; %Grav Const (m^3/(kg.s^2))
    ans = G*m_j*(r_j-r_i)/((norm(r_j-r_i))^3);
end
%%
function ans = dX_dT(m_j, r_i, r_j, v_i)
    if m_j == 0
        ans = 0;
    else
    ans = v_i;
    end
end