%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Four-Body Problem 
       %(Sun, Mercury, Earth, Mars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% sun Mercury Earth Mars 
%unit => Kg km kg s
%sun  1
[x1, v1] = planetEphemeris(juliandate(1969,7,16),'solarsystem','sun','432t','km');%AU/day
m1 = 1988500*10^24; %kg
x1 = double(x1'); v1 = double(v1');

%earth 2
[x2, v2] = planetEphemeris(juliandate(1969,7,16),'solarsystem','earth','432t','km');%AU/day
m2 = 5.97*10^24;
x2 = double(x2'); v2 = double(v2');

%mars 3 
[x3, v3] = planetEphemeris(juliandate(1969,7,16),'solarsystem','mars','432t','km');%AU/day
m3 = 0.642*10^24;
x3 = double(x3'); v3= double(v3');

%jupitor4
[x4, v4] = planetEphemeris(juliandate(1969,7,16),'solarsystem','Mercury','432t','km');%AU/day
m4 = 0.330*10^24;
x4 = double(x4'); v4= double(v4');


% constants
G = 6.67430*10^(-20); 
ti = 0.0000; %initial time
T = 1400*86400; %seconds %2 *orbital period of MARS
N = 3; %number of body 
S = 1400; % WE ARE VARYING THE STEP, NOT DT
dt = T/S; % timestep
% combine
u = [ x1; v1; x2; v2; x3; v3; x4; v4];
u0 =u;
%1:3 x1
%7:9 x2
%13:15 x3
%19:21 x4

%function
du = @(u) [...
    u(4);u(5);u(6);...
    G*((m2*(u(7:9)-u(1:3)))/(norm(u(7:9)-u(1:3)))^3 +...
    (m3*(u(13:15)-u(1:3)))/(norm(u(13:15)-u(1:3)))^3+...
    (m4*(u(19:21)-u(1:3)))/(norm(u(19:21)-u(1:3)))^3);...
    u(10);u(11);u(12);...
    G*((m1*(u(1:3)-u(7:9)))/(norm(u(1:3)-u(7:9)))^3 +...
    (m3*(u(13:15)-u(7:9)))/(norm(u(13:15)-u(7:9)))^3+...
    (m4*(u(19:21)-u(7:9)))/(norm(u(19:21)-u(7:9)))^3);...
    u(16);u(17);u(18);...
    G*((m1*(u(1:3)-u(13:15)))/(norm(u(1:3)-u(13:15)))^3 + ...
    (m2*(u(7:9)-u(13:15)))/(norm(u(7:9)-u(13:15)))^3+...
    (m4*(u(19:21)-u(13:15)))/(norm(u(19:21)-u(13:15)))^3);...
    u(22);u(23);u(24);...
    G*((m1*(u(1:3)-u(19:21)))/(norm(u(1:3)-u(19:21)))^3 + ...
    (m2*(u(7:9)-u(19:21)))/(norm(u(7:9)-u(19:21)))^3+...
    (m3*(u(13:15)-u(19:21)))/(norm(u(13:15)-u(19:21)))^3);...
    ];

%contruct approximation method
%forward euler
u_fe_keep = zeros(24,S);
u_fe_k = u;

%heun
u_h_keep = zeros(24,S);
u_h_k = u;

%RK4 
u_RK4_keep = zeros(24,S);
u_RK4_k = u;

%KDK method
x_KDK = [ x1; x2; x3; x4];
u_KDK = [ v1; v2; v3; v4];

%1:3 x1
%4:6 x2
%7:9 x3
%10:12 x4
du_KDK =@(X) [ ...
    G*((m2*(X(4:6)-X(1:3)))/(norm(X(4:6)-X(1:3)))^3 +...
    (m3*(X(7:9)-X(1:3)))/(norm(X(7:9)-X(1:3)))^3+...
    (m4*(X(10:12)-X(1:3)))/(norm(X(10:12)-X(1:3)))^3);...
    G*((m1*(X(1:3)-X(4:6)))/(norm(X(1:3)-X(4:6)))^3 +...
    (m3*(X(7:9)-X(4:6)))/(norm(X(7:9)-X(4:6)))^3+...
    (m4*(X(10:12)-X(4:6)))/(norm(X(10:12)-X(4:6)))^3);
    G*((m1*(X(1:3)-X(7:9)))/(norm(X(1:3)-X(7:9)))^3 + ...
    (m2*(X(4:6)-X(7:9)))/(norm(X(4:6)-X(7:9)))^3+...
    (m4*(X(10:12)-X(7:9)))/(norm(X(10:12)-X(7:9)))^3);...
    G*((m1*(X(1:3)-X(10:12)))/(norm(X(1:3)-X(10:12)))^3 + ...
    (m2*(X(4:6)-X(10:12)))/(norm(X(4:6)-X(10:12)))^3+...
    (m3*(X(7:9)-X(10:12)))/(norm(X(7:9)-X(10:12)))^3);...
    ];
vi_KDK = u_KDK;
xi_KDK = x_KDK;
u_KDK_keep = zeros(12,T/dt);

%Approximation method loops
tk = ti;
tic % timing the run of each method 

for i = 1 :S
    %KDK method 
    x_KDK_half = vi_KDK +(du_KDK(xi_KDK)*(dt/2));
    xi1= xi_KDK +  x_KDK_half*dt;
    vi1 = x_KDK_half + du_KDK(xi1)*dt/2;
    xi_KDK = xi1;
    vi_KDK = vi1;
    u_KDK_keep(:,i) = xi1; 
end 
t_KDK=toc;
tic
for i = 1:S 
    %forward euler
    u_fe_kp1 = u_fe_k +dt*du(u_fe_k);
    u_fe_k = u_fe_kp1;
    u_fe_keep(:,i) = u_fe_kp1;
    
end 
t_fe = toc;
tic
for i = 1 :S
    %heun's method
    u_h_kp1 = (u_h_k + 1/2*dt*(du(u_h_k)+ du(u_h_k + dt*du(u_h_k))));
    u_h_k = u_h_kp1;
    u_h_keep(:,i) = u_h_kp1;
    
end
t_h = toc;
tic

for i = 1 :S
    %RK4 
    y1 = du(u_RK4_k);
    y2 = du(u_RK4_k + dt*y1/2);
    y3 = du(u_RK4_k + dt*y2/2);
    y4 = du(u_RK4_k + dt*y3);
    
    u_RK4_kp1 = (u_RK4_k +1/6*dt*(y1+2*y2+2*y3+y4));
    u_RK4_k = u_RK4_kp1;
    u_RK4_keep(:,i) = u_RK4_kp1;
   
    tk=tk+dt;
end
 t_rk4 = toc; 



%% Obtain actual orbit location
steps = T/86400; %Total number of points
mtrue = zeros(steps,3);
etrue = zeros(steps,3);
Mtrue = zeros(steps,3);
for i = 1:steps
    [xM_t, vM_t] = planetEphemeris(juliandate(1969,7,16+(i)),'solarsystem','mars','432t','km');
    mtrue(i,:) = xM_t;
    [xE_t, vE_t] = planetEphemeris(juliandate(1969,7,16+(i)),'solarsystem','earth','432t','km');
    etrue(i,:) = xE_t;
     [xJ_t, vJ_t] = planetEphemeris(juliandate(1969,7,16+(i)),'solarsystem','Mercury','432t','km');
    Mtrue(i,:) = xJ_t;
end
sizextrue = size(mtrue); 


%% Plot

figure(1)

%fe method
subplot(2,2,1)
scatter3(u_fe_keep(1,:),u_fe_keep(2,:),u_fe_keep(3,:),5,'filled', 'r'), hold on
plot3(u_fe_keep(19,:),u_fe_keep(20,:),u_fe_keep(21,:),'linewidth',2)
plot3(u_fe_keep(7,:),u_fe_keep(8,:),u_fe_keep(9,:),'linewidth',2)
plot3(u_fe_keep(13,:),u_fe_keep(14,:),u_fe_keep(15,:),'linewidth',2)
% plot3(mtrue(:,1),mtrue(:,2),mtrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(etrue(:,1),etrue(:,2),etrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(Mtrue(:,1),Mtrue(:,2),Mtrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
%legend('sun','earth','true', 'mars')
title('\textbf{FE method}','interpreter', 'latex', 'fontsize', 18)
xlabel('x (km)','interpreter', 'latex', 'fontsize', 18)
ylabel('y (km)','interpreter', 'latex', 'fontsize', 18)
zlabel('z (km)','interpreter', 'latex', 'fontsize', 18)
legend('Sun','Mercury', 'Earth','Mars','interpreter','latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

% HEUN method
subplot(2,2,2) 
scatter3(u_fe_keep(1,:),u_fe_keep(2,:),u_fe_keep(3,:),5,'filled', 'r'), hold on
plot3(u_h_keep(7,:),u_h_keep(8,:),u_h_keep(9,:),'linewidth',2)
plot3(u_h_keep(13,:),u_h_keep(14,:),u_h_keep(15,:),'linewidth',2)
plot3(u_h_keep(19,:),u_h_keep(20,:),u_h_keep(21,:),'linewidth',2)
% plot3(mtrue(:,1),mtrue(:,2),mtrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(etrue(:,1),etrue(:,2),etrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(Mtrue(:,1),Mtrue(:,2),Mtrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
title("\textbf{Heun's Method}",'interpreter', 'latex', 'fontsize', 18)
xlabel('x (km)','interpreter', 'latex', 'fontsize', 18)
ylabel('y (km)','interpreter', 'latex', 'fontsize', 18)
zlabel('z (km)','interpreter', 'latex', 'fontsize', 18)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

%RK4 method
subplot(2,2,4)
scatter3(u_fe_keep(1,:),u_fe_keep(2,:),u_fe_keep(3,:),5,'filled', 'r'), hold on
plot3(u_RK4_keep(7,:),u_RK4_keep(8,:),u_RK4_keep(9,:),'linewidth',2)
plot3(u_RK4_keep(13,:),u_RK4_keep(14,:),u_RK4_keep(15,:),'linewidth',2)
plot3(u_RK4_keep(19,:),u_RK4_keep(20,:),u_RK4_keep(21,:),'linewidth',2)
% plot3(mtrue(:,1),mtrue(:,2),mtrue(:,3),'r--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(etrue(:,1),etrue(:,2),etrue(:,3),'b--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(Mtrue(:,1),Mtrue(:,2),Mtrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
title('\textbf{RK4 Method}','interpreter', 'latex', 'fontsize', 18)
xlabel('x (km)','interpreter', 'latex', 'fontsize', 18)
ylabel('y (km)','interpreter', 'latex', 'fontsize', 18)
zlabel('z (km)','interpreter', 'latex', 'fontsize', 18)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

%KDK method
subplot(2,2,3) 
scatter3(u_fe_keep(1,:),u_fe_keep(2,:),u_fe_keep(3,:),5,'filled', 'r'), hold on
plot3(u_KDK_keep(4,:),u_KDK_keep(5,:),u_KDK_keep(6,:),'linewidth',2)
plot3(u_KDK_keep(7,:),u_KDK_keep(8,:),u_KDK_keep(9,:),'linewidth',2)
plot3(u_KDK_keep(10,:),u_KDK_keep(11,:),u_KDK_keep(12,:),'linewidth',2)
% plot3(mtrue(:,1),mtrue(:,2),mtrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(etrue(:,1),etrue(:,2),etrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
% plot3(Mtrue(:,1),Mtrue(:,2),Mtrue(:,3),'--','LineWidth',1,'Color',[0.2, 0.4470, 0.410])
title('\textbf{KDK Method}','interpreter', 'latex', 'fontsize', 18)
xlabel('x (km)','interpreter', 'latex', 'fontsize', 18)
ylabel('y (km)','interpreter', 'latex', 'fontsize', 18)
zlabel('z (km)','interpreter', 'latex', 'fontsize', 18)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )


%% calculate radius at each point
E_true = double(etrue');
m_true = double(mtrue');
M_true = double(Mtrue');
r_E_true = zeros(steps,1);
r_m_true = zeros(steps,1);
r_M_true = zeros(steps,1);
for i = 1: steps
    r_E_true(i) = vecnorm(E_true(:,i));
    r_m_true(i) = vecnorm(m_true(:,i));
    r_M_true(i) = vecnorm(M_true(:,i));
end

%RK4
r_E_RK4 = zeros(S,1);
r_m_RK4 = zeros(S,1);
r_M_RK4 = zeros(S,1);
for i = 1: S
    r_E_RK4(i) = vecnorm(u_RK4_keep(7:9,i));
    r_m_RK4(i) = vecnorm(u_RK4_keep(13:15,i));
    r_M_RK4(i) = vecnorm(u_RK4_keep(19:21,i));
end

%forward euler
r_E_fe = zeros(S,1);
r_m_fe = zeros(S,1);
r_M_fe = zeros(S,1);
for i = 1: S
    r_E_fe(i) = vecnorm(u_fe_keep(7:9,i));
    r_m_fe(i) = vecnorm(u_fe_keep(13:15,i));
    r_M_fe(i) = vecnorm(u_fe_keep(19:21,i));
end

%heuns
r_E_h = zeros(S,1);
r_m_h = zeros(S,1);
r_M_h = zeros(S,1);
for i = 1: S
    r_E_h(i) = vecnorm(u_h_keep(7:9,i));
    r_m_h(i) = vecnorm(u_h_keep(13:15,i));
    r_M_h(i) = vecnorm(u_h_keep(19:21,i));
end

%KDK
r_E_KDK = zeros(S,1);
r_m_KDK = zeros(S,1);
r_M_KDK = zeros(S,1);

for i = 1: S
    r_E_KDK(i) = vecnorm(u_KDK_keep(4:6,i));
    r_m_KDK(i) = vecnorm(u_KDK_keep(7:9,i));
    r_M_KDK(i) = vecnorm(u_KDK_keep(10:12,i));
end

%% Plot Percentage difference
time = [ dt :dt : T];
figure(46)
plot(time, ((r_M_h-r_M_true)./r_M_true*100),'linewidth',2), hold on
plot(time,((r_M_KDK-r_M_true)./r_M_true*100),'linewidth',2)
plot(time, ((r_M_RK4-r_M_true)./r_M_true*100),'linewidth',2,'Color',[0.2, 0.4470, 0.410] )
plot(time,(r_M_true-r_M_true),'linewidth',2,'color',[0.9290 0.6940 0.1250])
legend("Heun's Method",'KDK','RK4','True Orbit','interpreter', 'latex', 'fontsize', 16)
xlabel('Time (s)','interpreter', 'latex', 'fontsize', 18)
ylabel('Percent difference (\%)','interpreter', 'latex', 'fontsize', 18)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 18 )
xlim([0 12*10^7])
grid on
mercury_norm = [...
    norm((r_M_h-r_M_true))/norm(r_M_true)*100;...
    norm(((r_M_KDK-r_M_true)))/norm(r_M_true)*100;...
    norm((r_M_RK4-r_M_true))/norm(r_M_true)*100]