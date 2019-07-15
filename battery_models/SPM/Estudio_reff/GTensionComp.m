function GTensionComp(X1, Y1, Y2, Y3, Y4, Y5)%, Y6, Y7)
%CREATEFIGURE(X1, Y1, X2, Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 07-Jun-2018 05:18:16

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(X1,Y1,'DisplayName','r_{eff} = 1','Color',[1 0 0]);

% Create plot
plot(X1,Y2,'DisplayName','r_{eff} = 1e-3','Color',[0 1 0]);

% Create plot
plot(X1,Y3,'DisplayName','r_{eff} = 1e-6','Color',[0 0 1]);

% Create plot
plot(X1,Y4,'DisplayName','r_{eff} = 1e-9','Color',[1 1 0]);

% Create plot
plot(X1,Y5,'DisplayName','r_{eff} = 1e18','Color',[1 0 1]);
% 
% % Create plot
% plot(X1,Y6,'DisplayName','r_{eff} = 1e6','Color',[0 1 1]);
% 
% % Create plot
% plot(X1,Y7,'DisplayName','r_{eff} = 1e9','Color',[0 0 0]);

% Create xlabel
xlabel('Tiempo, t [s]');

% Create title
title('Comparaci�n para distintos r_{eff}');

% Create ylabel
ylabel('Tensi�n descarga SPM, V [V]');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14);
% Create legend
legend(axes1,'show');
