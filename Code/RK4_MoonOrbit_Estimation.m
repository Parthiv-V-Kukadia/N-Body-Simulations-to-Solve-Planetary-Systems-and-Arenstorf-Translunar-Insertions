    %% Initialize Values

clc
clear all
close all
%Earth %Moon %Other Satellites/ Capsules

N = 2; %number of bodies (<= 5)

Time = 660; %Total time elapsed (hours)
steps = 100000; %Total number of points

%Change up to the Nth Initial Body Position

[xE_i, vE_i] = planetEphemeris(juliandate(1969,7,16),'EarthMoon','Earth');
[xM_i, vM_i] = planetEphemeris(juliandate(1969,7,16),'EarthMoon','Moon');

%Find Proper Rotation

z_rot_vec = [xM_i(1) xM_i(2) 0];
y_rot_vec = [xM_i(1) 0 xM_i(3)];
x_rot_vec = [0 xM_i(2) xM_i(3)];
ang_x = acosd(dot(xM_i, x_rot_vec)/norm(xM_i)/norm(x_rot_vec));
ang_y = acosd(dot(xM_i, y_rot_vec)/norm(xM_i)/norm(y_rot_vec));
ang_z = acosd(dot(xM_i, z_rot_vec)/norm(xM_i)/norm(z_rot_vec));

syms x y z real
Rz = [cosd(z) -sind(z) 0; sind(z) cosd(z) 0; 0 0 1];
Ry = [cosd(y) 0 sind(y); 0 1 0; -sind(y) 0 cosd(y)];
Rx = [1 0 0; 0 cosd(x) -sind(x); 0 sind(x) cosd(x)];
Rot = subs(Rx*Rz,[x,z],[ang_x, ang_z]);

xE_i = 1000*xE_i;
vE_i = 1000*vE_i;
xM_i = 1000*xM_i;
vM_i = 1000*vM_i;

x_b1 = xE_i + xE_i/norm(xE_i);%*sqrt(1000000); %Initial 1st Body Position (m)
x_b2 = [3e7 3e7]; %Initial 2nd Body Position (m)
x_b3 = [0 0]; %Initial 3rd Body Position (m)

%Change up to the Nth Initial Body Velocity

v_b1 = [0 0 0]; %Initial 1st Body Velocity (m/s)
v_b2 = [50 0]; %Initial 2nd Body Velocity (m/s)
v_b3 = [0 0]; %Initial 3rd Body Velocity (m/s)

%Change up to the Nth Initial Body Mass
m_E = 5.97237e24; %Earth Mass (kg)
m_M = 7.34767309e22; %Moon Mass (kg)

m_b1 = 20000; %Initial 1st Body Mass (kg)
m_b2 = 2.972e7; %Initial 2nd Body Mass (kg)
m_b3 = 0; %Initial 3rd Body Mass (kg)
    
%%
x = zeros(steps, N*3);
v = zeros(steps, N*3);

xpos_m = [0 0 0];
vpos_m = [0 0 0];

xpos_e = [0 0 0];
vpos_e = [0 0 0];


if N == 2
    mass = [m_E, m_M];
    x(1,:) = [xE_i, xM_i];
    v(1,:) = [vE_i, vM_i];
elseif N == 3
    mass = [m_E, m_M, m_b1];
    x(1,:) = [xE_i, xM_i, x_b1];
    v(1,:) = [vE_i, vM_i, v_b1];
elseif N ==4
    mass = [m_E, m_M, m_b1, m_b2];
    x(1,:) = [xE_i, xM_i, x_b1, x_b2];
    v(1,:) = [vE_i, vM_i, v_b1, v_b2];
else
    mass = [m_E, m_M, m_b1, m_b2, m_b3];
    x(1,:) = [xE_i, xM_i, x_b1, x_b2, x_b3];
    v(1,:) = [vE_i, vM_i, v_b1, v_b2, v_b3];
end


t = 0;
dt = Time*3600/steps; %time step (s)


for k = 1:steps
    t = t + dt;
    for i = 1:N
        x_i = x(k, ((i-1)*3+1):i*3); %Get the position in m
        v_i = v(k, ((i-1)*3+1):i*3); %Get the vel in m/s
        
        x_ip1 = x_i;
        v_ip1 = v_i;
        const = 0;
        for j = 1:N
            if j ~= i
                x_j = x(k,((j-1)*3+1):j*3);
                k1x = dx_dt(mass(1,j), x_i, x_j, v_i);
                k1v = dv_dt(mass(1,j), x_i, x_j);
                k2x = dx_dt(mass(1,j), x_i+dt*k1x/2, x_j, v_i+dt*k1v/2);
                k2v = dv_dt(mass(1,j), x_i + dt*k1x/2, x_j);
                k3x = dx_dt(mass(1,j), x_i+dt*k2x/2, x_j, v_i+dt*k2v/2);
                k3v = dv_dt(mass(1,j), x_i + dt*k2x/2, x_j);
                k4x = dx_dt(mass(1,j), x_i+dt*k3x, x_j, v_i+dt*k3v);
                k4v = dv_dt(mass(1,j), x_i + dt*k3x, x_j);
                
                x_ip1 = x_ip1 + dt/6*(k1x+2*k2x+2*k3x+k4x);
                v_ip1 = v_ip1 + dt/6*(k1v+2*k2v+2*k3v+k4v);
            end
        end
        x(k+1, ((i-1)*3+1):i*3) = x_ip1;
        v(k+1, ((i-1)*3+1):i*3) = v_ip1;
        if k > 1
            if N ==2
                if x_i(1) < 10e2
                    xpos_m = [x_i(1) x_i(2) x_i(3)];
                    vpos_m = [v(k-1,4) v(k-1,5) v(k-1,6)];
                    xpos_e = [x(k-1,1) x(k-1,2) x(k-1,3)];
                    vpos_e = [v(k-1,1) v(k-1,2) v(k-1,3)];
                end
            end
        end
    end
end

%% 3D View
figure(1)
for i=1:N
    if i == 1
        plot3(x(:,(i-1)*3+1),x(:,(i-1)*3+2),x(:,(i-1)*3+3),'LineWidth',3)
        hold on
    end
    if i ~=3
        if i ~=1
            plot3(x(:,(i-1)*3+1),x(:,(i-1)*3+2),x(:,(i-1)*3+3),'--','LineWidth',3)
            hold on
        end
    else
        plot3(x(:,(i-1)*3+1),x(:,(i-1)*3+2),x(:,(i-1)*3+3),'k','--','LineWidth',3)
        hold on
    end
    if i==1
        set(gcf, 'PaperUnits', 'centimeters')
        set(gcf, 'PaperSize', [18 18])
        set(gcf, 'Units', 'centimeters' )
        set(gcf, 'Position', [8 2 18 18])
        set(gcf, 'PaperPosition', [0 0 18 18])
    end
end
xlabel('Distance (Km)')
ylabel('Distance (Km)')
zlabel('Distance (Km)')
legend('Earth','Moon(approx)','Moon(true)')
%% True Data
xtrue = zeros(steps/100,3);
for i = 1:steps/100
    [xM_t, vM_t] = planetEphemeris(juliandate(1969,7,27.922*100/steps*(i-1)+1),'EarthMoon','Moon');
    xtrue(i,:) = xM_t*1000;
end
%%
plot3(xtrue(:,1),xtrue(:,2),xtrue(:,3),'-','LineWidth',3,'Color',[0.2, 0.4470, 0.410])
%% Planar

xE_i = xE_i*Rot;
vE_i = vE_i*Rot;
xM_i = xM_i*Rot;
vM_i = vM_i*Rot;


if N == 2
    mass = [m_E, m_M];
    x(1,:) = [xE_i, xM_i];
    v(1,:) = [vE_i, vM_i];
elseif N == 3
    mass = [m_E, m_M, m_b1];
    x(1,:) = [xE_i, xM_i, x_b1];
    v(1,:) = [vE_i, vM_i, v_b1];
elseif N ==4
    mass = [m_E, m_M, m_b1, m_b2];
    x(1,:) = [xE_i, xM_i, x_b1, x_b2];
    v(1,:) = [vE_i, vM_i, v_b1, v_b2];
else
    mass = [m_E, m_M, m_b1, m_b2, m_b3];
    x(1,:) = [xE_i, xM_i, x_b1, x_b2, x_b3];
    v(1,:) = [vE_i, vM_i, v_b1, v_b2, v_b3];
end

t = 0;
dt = Time*3600/steps; %time step (s)


for k = 1:steps
    t = t + dt;
    for i = 1:N
        x_i = x(k, ((i-1)*3+1):i*3); %Get the position in m
        v_i = v(k, ((i-1)*3+1):i*3); %Get the vel in m/s
        
        x_ip1 = x_i;
        v_ip1 = v_i;
        const = 0;
        for j = 1:N
            if j ~= i
                x_j = x(k,((j-1)*3+1):j*3);
                k1x = dx_dt(mass(1,j), x_i, x_j, v_i);
                k1v = dv_dt(mass(1,j), x_i, x_j);
                k2x = dx_dt(mass(1,j), x_i+dt*k1x/2, x_j, v_i+dt*k1v/2);
                k2v = dv_dt(mass(1,j), x_i + dt*k1x/2, x_j);
                k3x = dx_dt(mass(1,j), x_i+dt*k2x/2, x_j, v_i+dt*k2v/2);
                k3v = dv_dt(mass(1,j), x_i + dt*k2x/2, x_j);
                k4x = dx_dt(mass(1,j), x_i+dt*k3x, x_j, v_i+dt*k3v);
                k4v = dv_dt(mass(1,j), x_i + dt*k3x, x_j);
                
                x_ip1 = x_ip1 + dt/6*(k1x+2*k2x+2*k3x+k4x);
                v_ip1 = v_ip1 + dt/6*(k1v+2*k2v+2*k3v+k4v);
            end
        end
        x(k+1, ((i-1)*3+1):i*3) = x_ip1;
        v(k+1, ((i-1)*3+1):i*3) = v_ip1;
    end
end

%% Planar View
figure(2)
for i=1:N
    if i == 1
        plot3(x(:,(i-1)*3+1),x(:,(i-1)*3+2),x(:,(i-1)*3+3),'LineWidth',3)
        hold on
    end
    if i ~=3
        if i ~=1
            plot3(x(:,(i-1)*3+1),x(:,(i-1)*3+2),x(:,(i-1)*3+3),'--','LineWidth',3)
            hold on
        end
    else
        plot3(x(:,(i-1)*3+1),x(:,(i-1)*3+2),x(:,(i-1)*3+3),'k','--','LineWidth',3)
        hold on
    end
    if i==1
        set(gcf, 'PaperUnits', 'centimeters')
        set(gcf, 'PaperSize', [18 18])
        set(gcf, 'Units', 'centimeters' )
        set(gcf, 'Position', [8 2 20 20])
        set(gcf, 'PaperPosition', [0 0 20 20])
        set(gca,'interpreter','latex','fontsize',16)
        view(2)
    end
end
xnew = xtrue*Rot;
plot3(xnew(:,1),xnew(:,2),xnew(:,3),'-','LineWidth',3,'Color',[0.2, 0.4470, 0.410])

xlabel('Distance (km)')
ylabel('Distance (km)')
legend('Earth','Moon (RK4)','Moon (True)')
%%
function ans = dv_dt(m_j, r_i, r_j)
    G = 6.67430e-11; %Grav Const (m^3/(kg.s^2))
    ans = G*m_j*(r_j-r_i)/((norm(r_j-r_i))^3);
end
%%
function ans = dx_dt(m_j, r_i, r_j, v_i)
    if m_j == 0
        ans = 0;
    else
    ans = v_i;
    end
end