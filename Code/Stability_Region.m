%Stability plots for 4 different numerical methods:
%Heun's method, RK4, Forward Euler, and Leapfrog/Midpoint 

%Clear all previous variables, and close all previous windows
clear all; close all; clc;

% Specifying the x and y range that will be used for plotting the graphs
xl = -4; xr = 4;
x_axis = [xl,xr];
yl = -4; yr = 4;
y_axis = [yl,yr];

% Constructing a mesh for the graph
x_lin = linspace(xl,xr,401);
y_lin = linspace(yl,yr,401);
[x,y] = meshgrid(x_lin,y_lin);

% Calculating z for z-transform domain
z = x + y*i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stability calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Heun's method
Heun = 1 + z + 0.5*z.^2;
Heun = abs(Heun);

%RK4 
RK4 = 1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4;
RK4 = abs(RK4);

%Forward Euler
t = linspace(0,2*pi,200); %Theta
l_d_t = exp(t*i)-1; %lamda * delta t
x_R = real(l_d_t); %Real part of lambda * delta t
y_I = imag(l_d_t); %Imaginary part of lambda * delta t

%Leapfrog
%Stability region for Leapfrog is only on the imaginary axis, and not on
%the real axis, therefore, it is just a vertical line within the limits of
%-1 and 1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Stability Regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

%Plotting Heun's stability
contour(x,y,Heun,[1 1],'color','[0.8500, 0.3250, 0.0980].','linewidth',3);
hold on;

%Plotting RK4 Stability
contour(x,y,RK4,[1 1],'b--','linewidth',3);
hold on;

%Plotting Forward Euler Stability
plot(x_R,y_I,'m-','linewidth',3)
hold on;

%Plotting Leapfrog Stability
line([0,0],[-1,1],'linestyle','-','color','r','linewidth',4);

%Setting the graph up
axis([x_axis,y_axis]);
grid on;

%Graph Labels
xlabel('Real($\lambda_{i} \Delta t$)','interpreter','latex');
ylabel('Imaginary($\lambda_{i} \Delta t$)','interpreter','latex');
legend({'Heuns','RK4','Forward Euler','Leapfrog'},'interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',18)

%Scaling
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [22 22])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 22 22])
set(gcf, 'PaperPosition', [0 0 22 22])
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,'Stability.jpg')