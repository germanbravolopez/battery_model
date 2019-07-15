% Modelo de batería de Li-ion basado en el modelo SPM que utiliza la
% función "pdepe" para integrar en el espacio y el tiempo la ecuación de la
% concentración de litio en la partícula sólida.

% --------------------------------------------------------------
% --------------------------------------------------------------
% Parámetros de simulación
% --------------------------------------------------------------

i_salida = [5.26852693329723];%[5.26823479036215]; %4.72772522522523 23.517; % [A] Corriente de salida de la celda
r = linspace(0,1e-4,60); % [cm] La distancia completa del electrodo
t = linspace(0,9000,6000); % [s] Tiempo de simulación 3930

% --------------------------------------------------------------
% Parámetros de la tabla FreedomCar
% --------------------------------------------------------------
tic
% Design specifications (geometry and volume fractions)
delta_p = 36.4e-4; delta_n = 50e-4; % [cm] Thickness
L_p = 36.4e-4; L_n = 50e-4; % [cm] Thickness
R_p = 1e-4; % [cm] Particle radius
epsilon_sp = 0.5; epsilon_sn = 0.58; % Active material volume fraction
A = 10452; % [cm^2] Electrode plate area
global I
I = i_salida/A; % [A*cm^-2] Corriente por unidad de área

% Solid and electrolyte phase Li-ion concentration (0% and 100% SOC)
global beta_n100 beta_n0 beta_p100 beta_p0 c_smaxp c_smaxn
c_smaxp = 23.9e-3; c_smaxn = 16.1e-3; % [mol*cm^-3] Maximum solid phase concentration
beta_n100 = 0.676; beta_n0 = 0.126; % Negative electrode stoichiometry 
beta_p100 = 0.442; beta_p0 = 0.936; % Positive electrode stoichiometry
c_e = 1.2e-3; % [mol*cm^-3] Initial electrolyte phase concentration

% Kinetic and transport properties
i_0p = 2.6e-3; i_0n = 3.6e-3; % [A*cm^-2] Exchange current density
D_sp = 3.7e-12; D_sn = 2e-12; % [cm^2*s^-1] Solid phase Li diffusion coeficient
D_e = 2.6e-6; % [cm^2*s^-1] Electrolyte phase Li-ion diffusion coeficient

R_collector = 1.9e-3; % [Ohms] Current collector contact resistant
R_f = 0; % [Ohm*cm^-1] Solid-electrolyte interphase layer film resistance

% Constantes
F = 96478; % [C*mol^-1] Faraday's constant
R = 8.314; % [J*mol^-1] Universal gas constant 
T = 298; % [K] Temperature
r_eff = 1; % [mol*cm^-2*s^-1*(mol*cm^-1*s)^-1.5] Reaction rate constant

% --------------------------------------------------------------
% Ecuaciones del modelo
% --------------------------------------------------------------

% Equation of the specific area
    a_p = 3*epsilon_sp/R_p;
    a_n = 3*epsilon_sn/R_p;

% Molar flux equation
    j_np = -(I)./(F*a_p*L_p);
    j_nn = (I)./(F*a_n*L_n);

% Li-ion concentration in the solid particles
%     global cont cont1 cont2
%     cont = 1e9; cont1 = 0; cont2 = 0;
%     c_sp_prueba = concentration_p(r,t);
%     I = i_salida/A;
%     c_sn_prueba = concentration_n(r,t);
%     I = i_salida/A;
%     cont = min(cont1, cont2);
    
%     cont1 = 0; 
    c_sp = concentration_p(r,t);
    I = i_salida/A;
    
%     cont2 = 0;
    c_sn = concentration_n(r,t);
    I = i_salida/A;
    

    
% Surface concentration of Li-ion in the solid particle
    c_ssp = c_sp(:,end);
    c_ssn = c_sn(:,end);
    
% Stoichiometry for the potential
    beta_p = c_ssp./c_smaxp;
    beta_n = c_ssn./c_smaxn;
    
% Open circuit potential: p->(+)  n->(-)
    U_p = 85.681.*beta_p.^6-357.7.*beta_p.^5+613.89.*beta_p.^4-555.65.*beta_p.^3+281.06.*beta_p.^2-76.648.*beta_p-0.30987.*exp(5.657.*beta_p.^115)+13.1983;
    U_n = 8.00229+5.0647.*beta_n-12.578.*beta_n.^(1/2)-8.6322e-4.*beta_n.^(-1)+2.1765e-5.*beta_n.^(3/2)-0.46016.*exp(15.*(0.06-beta_n))-0.55364.*exp(-2.4326.*(beta_n-0.92));

% Potential equations
    theta_s0p = (2.*R.*T./F).*asinh(-I./(2.*a_p.*L_p.*r_eff.*sqrt(c_e.*c_ssp.*(c_smaxp-c_ssp)))) + U_p + (R_f.*I)./(a_p.*L_p);
    theta_s0n = (2.*R.*T./F).*asinh(I./(2.*a_n.*L_n.*r_eff.*sqrt(c_e.*c_ssn.*(c_smaxn-c_ssn)))) + U_n + (R_f.*I)./(a_n.*L_n);
    
% Output voltage
    V = theta_s0p - theta_s0n - I.*R_collector;

    figure(4)
    plot(t,V)
    title('Tensión de salida de la batería')
    xlabel('Tiempo, t [s]')
    ylabel('Tensión, V [V]')

% --------------------------------------------------------------
% Bring C-rate off
    [x_p y_p] = find(c_sp == min(min((c_sp))));
    [x_n y_n] = find(c_sn == max(max((c_sn))));
    Q_batt = abs(i_salida*t(x_p)/3600); % [Ah] Carga de la batería
    Q_batt1 = abs(i_salida*t(x_n)/3600);
    C_rate = abs(i_salida/Q_batt); % Rate of discharge
    

%     x = min(x_p,x_n)


tSim = toc % [s] Tiempo de que emplea el programa en ejecutarse