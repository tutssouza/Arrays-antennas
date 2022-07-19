%ARTHUR SILVA SOUZA

clear all; 
close all; 
clc

%% Array factor

% Valores de entrada

lambda = 1; % Valor arbitrario para essa aplicacao
beta = -pi/2;
d = lambda/4;
k = 2*pi/lambda;
theta = 0:pi/100:2*pi;

AF = abs((cos((1/2)*(k*d*cos(theta) + beta))));
AF_dB = 20*log10(AF);

figure (1);

subplot(2,2,1);
plot (theta,AF,'LineWidth',3)
title('Array factor ','FontSize',14);
xlabel('\phi','fontweight','bold','FontSize',14);
xlim([0 2*pi]);
ylabel('Amplitude','fontweight','bold','FontSize',14);

subplot(2,2,2);
plot (theta,AF_dB,'LineWidth',3)
title('Array factor (in dB)','FontSize',14);
xlabel('\phi','fontweight','bold','FontSize',14);
xlim([0 2*pi]);
ylabel('Amplitude','fontweight','bold','FontSize',14);

% plotando de forma polar

rmin = min(AF);
rmax = max(AF);

subplot(2,2,3);
polarplot (AF, 'LineWidth',3)
rlim([rmin rmax]);
ax = gca;
ax.ThetaZeroLocation = 'top';
plotTitle = sprintf('Radiation Pattern (normalized)');
title(plotTitle)

% plotando de forma polar em dB
rmin_dB = min(AF_dB);
rmax_dB = max(AF_dB);

subplot(2,2,4);
polarplot (AF_dB, 'LineWidth',3);
rlim([rmin_dB rmax_dB]);
ax = gca;
ax.ThetaZeroLocation = 'top';
plotTitle = sprintf('Radiation Pattern (in dB)');
title(plotTitle);

%% Element 

% Valores de entrada

lambda = 1; % Valor arbitrario para essa aplicacao
beta = 0;
d = lambda/4;
k = 2*pi/lambda;
theta = 0:pi/100:2*pi;

% Definindo equacao

E = abs(2*cos(theta).*(cos((1/2)*(k*d*cos(theta) + beta))));
E_dB = 20*log10(E);

figure (2);

subplot(2,2,1);
plot (theta,E,'LineWidth',3)
title('Element','FontSize',14);
xlabel('\phi','fontweight','bold','FontSize',14);
xlim([0 2*pi]);
ylabel('Amplitude','fontweight','bold','FontSize',14);

subplot(2,2,2);
plot (theta,E_dB,'LineWidth',3)
title('Element (in dB)','FontSize',14);
xlabel('\phi','fontweight','bold','FontSize',14);
xlim([0 2*pi]);
ylabel('Amplitude','fontweight','bold','FontSize',14);

% plotando de forma polar

rmin = min(E);
rmax = max(E);

subplot(2,2,3);
polarplot (E, 'LineWidth',3)
rlim([rmin rmax]);
ax = gca;
ax.ThetaZeroLocation = 'top';
plotTitle = sprintf('Radiation Pattern (normalized)');
title(plotTitle)

% plotando de forma polar em dB
rmin_dB = min(E_dB);
rmax_dB = max(E_dB);

subplot(2,2,4);
polarplot (E_dB, 'LineWidth',3);
rlim([rmin_dB rmax_dB]);
ax = gca;
ax.ThetaZeroLocation = 'top';
plotTitle = sprintf('Radiation Pattern (in dB)');
title(plotTitle);

%% Equally spaced linear array

% Valores de entrada

N = 5;
phi = 0:pi/100:2*pi;

f = abs(sin(N*phi/2)./(N*sin(phi/2)));

figure (3);

subplot(2,2,4);
plot (phi,f,'LineWidth',3)
title('Array factor for an equally spaced,','FontSize',14);
xlabel('\phi','fontweight','bold','FontSize',14);
xlim([0 2*pi]);
ylabel('Amplitude','fontweight','bold','FontSize',14);

subplot(2,2,3);
polarplot (f, 'LineWidth',3)
rmin = min(f);
rmax = max(f);
rlim([rmin rmax]);
ax = gca;
ax.ThetaZeroLocation = 'top';
plotTitle = sprintf('Radiation Pattern (normalized)');
title(plotTitle)

% Plotando em 3D

xv = f.* cos(phi);  %  transformando para polar
yv = f.* sin(phi);  %  transformado para polar
  
angulo_de_revolucao = 0:0.1:2*pi; %  Angulo de revolucao

% Criando pontos em 3D

xf = repmat(xv',size(angulo_de_revolucao)); 
yf = yv' * cos(angulo_de_revolucao);
zf = yv' * sin(angulo_de_revolucao);

subplot(2,2,1);
mesh(xf,yf,zf);
plotTitle = sprintf('3D Radiation Pattern (normalized)');
title(plotTitle)

subplot(2,2,2);
mesh(zf,yf,xf);
view(90,0)
plotTitle = sprintf('3D Radiation Pattern (normalized)');
title(plotTitle)

%% Broadside Array

% Valores de entrada

N = 10;
lambda = 1;
beta = 0;
d = lambda;
k = 2*pi/lambda;
theta = 0:pi/100:2*pi;

phi_m = k*d*cos(theta)+beta;
f = abs(sin(N*phi_m/2)./(N*sin(phi_m/2)));

figure (4);

subplot(2,2,4);
plot (phi,f,'LineWidth',3)
title('Array factor for an equally spaced,','FontSize',14);
xlabel('\phi','fontweight','bold','FontSize',14);
xlim([0 2*pi]);
ylabel('Amplitude','fontweight','bold','FontSize',14);

subplot(2,2,3);
polarplot (f, 'LineWidth',3)
rmin = min(f);
rmax = max(f);
rlim([rmin rmax]);
ax = gca;
ax.ThetaZeroLocation = 'top';
plotTitle = sprintf('Radiation Pattern (normalized)');
title(plotTitle)

% Plotando em 3D

xv = f.* cos(phi);  %  transformando para polar
yv = f.* sin(phi);  %  transformado para polar
  
angulo_de_revolucao = 0:0.1:2*pi; %  Angulo de revolucao

% Criando pontos em 3D

xf = repmat(xv',size(angulo_de_revolucao)); 
yf = yv' * cos(angulo_de_revolucao);
zf = yv' * sin(angulo_de_revolucao);

subplot(2,2,1);
mesh(xf,yf,zf);
plotTitle = sprintf('3D Radiation Pattern (normalized)');
title(plotTitle)

subplot(2,2,2);
mesh(zf,yf,xf);
view(90,0)
plotTitle = sprintf('3D Radiation Pattern (normalized)');
title(plotTitle)