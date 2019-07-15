clear all
close all

V_spm = load('V_spm.txt','V','-ascii');
V_3p = load('V_3p.txt','V','-ascii');
t_spm = load('t_spm.txt','t','-ascii');
t_3p = load('t_3p.txt','t','-ascii');


figure(1)
plot(t_spm,V_spm,'r')
hold on
plot(t_3p,V_3p,'b')
title('Comparación de tensiones SPM y Tres-parámetros')
xlabel('Tiempo, t [s]')
ylabel('Tensión, V [V]')
legend('1C-Tensión salida SPM','1C-Tensión salida Tres-parámetros')