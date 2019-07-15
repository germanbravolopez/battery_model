clear all
close all

t_spm = load('t_spm.txt','t','-ascii');
V_spm = load('V_spm.txt','V','-ascii');
V_spm_r1 = load('V_spm_r1.txt','V','-ascii');
V_spm_r2 = load('V_spm_r2.txt','V','-ascii');
V_spm_r3 = load('V_spm_r3.txt','V','-ascii');
V_spm_r4 = load('V_spm_r4.txt','V','-ascii');
V_spm_r5 = load('V_spm_r5.txt','V','-ascii');
V_spm_r6 = load('V_spm_r6.txt','V','-ascii');
V_spm_r7 = load('V_spm_r7.txt','V','-ascii');


GTensionComp(t_spm, V_spm, V_spm_r1, V_spm_r2, V_spm_r3, V_spm_r7)
     
GTensionComp(t_spm, V_spm,V_spm_r4, V_spm_r5, V_spm_r6, V_spm_r7)