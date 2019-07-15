% Li-ion battery model based on the Three parameters model. Which is
% solved from the SPM with a volume average.
clear all
close all
clf
% --------------------------------------------------------------
% --------------------------------------------------------------
% Simulation parameters
% --------------------------------------------------------------
i_salida = 6.0104;        % [A] Cell's output current
t = linspace(0,4500,800); % [s] Simulation time

% --------------------------------------------------------------
% FreedomCar table's parameters 
% --------------------------------------------------------------
tic
% Design specifications (geometry and volume fractions)
delta_p = 36.4e-4; delta_n = 50e-4;  % [cm] Thickness
L_p = 36.4e-4; L_n = 50e-4;          % [cm] Thickness
R_p = 1e-4;                          % [cm] Particle radius
epsilon_sp = 0.5; epsilon_sn = 0.58; % Active material volume fraction
A = 10452;                           % [cm^2] Electrode plate area
I = i_salida/A;                      % [A*cm^2] Current per area unit

% Solid and electrolyte phase Li-ion concentration (0% and 100% SOC)
c_smaxp = 23.9e-3; c_smaxn = 16.1e-3; % [mol*cm^-3] Maximum solid phase concentration
beta_n100 = 0.676; beta_n0 = 0.126;   % Negative electrode stoichiometry 
beta_p100 = 0.442; beta_p0 = 0.936;   % Positive electrode stoichiometry
c_e = 1.2e-3; % [mol*cm^-3] Initial electrolyte phase concentration

% Kinetic and transport properties
i_0p = 2.6e-3; i_0n = 3.6e-3; % [A*cm^-2] Exchange current density
D_sp = 3.7e-12; D_sn = 2e-12; % [cm^2*s^-1] Solid phase Li diffusion coeficient
D_e = 2.6e-6; % [cm^2*s^-1] Electrolyte phase Li-ion diffusion coeficient

R_collector = 1.9e-3; % [Ohms] Current collector contact resistant
R_SEI = 0; % [Ohm*cm^-1] Solid-electrolyte interphase layer film resistance

% Constants
F = 96478; % [C*mol^-1] Faraday's constant
R = 8.314; % [J*mol^-1] Universal gas constant 
T = 298;   % [K] Temperature
r_eff = 1; % [mol*cm^-2*s^-1*(mol*cm^-1*s)^-1.5] Reaction rate constant

% --------------------------------------------------------------
% Model equations
% --------------------------------------------------------------

% Equation of the specific area
    a_p = 3*epsilon_sp/R_p;
    a_n = 3*epsilon_sn/R_p;

for i = 1:1:length(t)
% Molar flux equation
    j_np = -(I)/(F*a_p*L_p);
    j_nn = (I)/(F*a_n*L_n);
    
% Volume-average solid phase concentration
    c_sp(i) = -(3*j_np.*t(i))./(R_p) + beta_p100*c_smaxp; % Battery charged
    c_sn(i) = -(3*j_nn.*t(i))./(R_p) + beta_n100*c_smaxn; % at t=0
    
% Volume-averaged concentration flux
    q_sp = j_np*(-45/(2*R_p^2))/(30*D_sp/R_p^2); % Integral's solution 
    q_sn = j_nn*(-45/(2*R_p^2))/(30*D_sn/R_p^2); 
    
% Surface concentration of lithium in the particles
    c_ssp(i) = c_sp(i) + (8*R_p*q_sp)/(35) - (R_p.*j_np)./(35*D_sp);
    c_ssn(i) = c_sn(i) + (8*R_p*q_sn)/(35) - (R_p.*j_nn)./(35*D_sn);
    
% Stoichiometry for the potential
    beta_p(i) = c_ssp(i)./c_smaxp;
    if (beta_p(i) > beta_p0) || (beta_p(i) < beta_p100)
        I = 0;
    end
    beta_n(i) = c_ssn(i)./c_smaxn;
    if (beta_n(i) > beta_n100) || (beta_n(i) < beta_n0)
        I = 0;
    end
    
% Open circuit potential
    U_p(i) = 85.681.*beta_p(i).^6-357.7.*beta_p(i).^5+613.89.*beta_p(i).^4-555.65.*beta_p(i).^3+281.06.*beta_p(i).^2-76.648.*beta_p(i)-0.30987.*exp(5.657.*beta_p(i).^115)+13.1983;
    U_n(i) = 8.00229+5.0647.*beta_n(i)-12.578.*beta_n(i).^(1/2)-8.6322e-4.*beta_n(i).^(-1)+2.1765e-5.*beta_n(i).^(3/2)-0.46016.*exp(15.*(0.06-beta_n(i)))-0.55364.*exp(-2.4326.*(beta_n(i)-0.92));

% Potential equations
    theta_s0p(i) = (2.*R.*T./F).*asinh(-I./(2.*a_p.*L_p.*r_eff.*sqrt(c_e.*c_ssp(i).*(c_smaxp-c_ssp(i))))) + U_p(i) + (R_SEI.*I)./(a_p.*L_p);
    theta_s0n(i) = (2.*R.*T./F).*asinh(I./(2.*a_n.*L_n.*r_eff.*sqrt(c_e.*c_ssn(i).*(c_smaxn-c_ssn(i))))) + U_n(i) + (R_SEI.*I)./(a_n.*L_n);
    
% Output voltage
    V(i) = theta_s0p(i) - theta_s0n(i) - I.*R_collector;
end

figure(1)
plot(t,V)
title('Tensión de salida de la batería')
xlabel('Tiempo, t [s]')
ylabel('Tensión, V [V]')

figure(2)
subplot(2,1,1)
    plot(t,c_sp)
    title('Concentración electrodo^{+}')
    ylabel('c_{sp} [mol*cm^{-3}]')
    xlabel('Tiempo, t [s]')
subplot(2,1,2)
    plot(t,c_sn)
    title('Concentración electrodo^{-}')
    ylabel('c_{sn} [mol*cm^{-3}]')
    xlabel('Tiempo, t [s]')
    
figure(3)
subplot(2,1,1)
    plot(beta_p,U_p)
    title('Potencial^{+} circuito abierto')
    xlabel('Coeficiente estequiométrico, beta_p')
    ylabel('Potencial, U_{+} [V]')
subplot(2,1,2)
    plot(beta_n,U_n)
    title('Potencial^{-} circuito abierto')
    xlabel('Coeficiente estequiométrico, beta_n')
    ylabel('Potencial, U_{-} [V]')
    
% --------------------------------------------------------------
% Calculate C-rate
    [x_p, y_p] = find(c_sp == max(c_sp));
    Q_batt = abs(i_salida*t(y_p)/3600);  % [Ah] Battery capacity
    C_rate = abs(i_salida/Q_batt)        % Rate of discharge

    
tSim = toc  % [s] Time the program takes to end

% savdir = 'C:\Users\Usuario\Desktop\German\UNIVERSIDAD\4o\2o_cuatrimestre\TFG\Modelos_batt\BATT_MODEL';
% save(fullfile(savdir,'V_3p.txt'),'V','-ascii');
% save(fullfile(savdir,'t_3p.txt'),'t','-ascii');