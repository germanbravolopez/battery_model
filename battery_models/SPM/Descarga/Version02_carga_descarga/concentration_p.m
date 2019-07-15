function [c_sp] = concentration_p(r,t)

m = 2;
x=r;

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);
c_sp = u;

% A surface plot is often a good way to study a solution.
figure(1)
surf(x,t,u) 
title('Concentration electrode^{+}.')
xlabel('Distance, r [cm]')
ylabel('Time, t [s]')
zlabel('Concentration^+, c_{s} [mol*m^{-3}]')

% A solution profile can also be illuminating.
figure(3)
subplot(2,1,1);
plot(t,u(:,end))
title('Concentration electrode^{+} at R_{p}')
xlabel('Time, t [s]')
ylabel('c_{s}(R_{p},t) [mol*cm^{-3}]')
% --------------------------------------------------------------
function [c,f,s] = pdefun(x,t,u,DuDx)
D_s = 3.7e-12;

c = 1/D_s;
f = DuDx;
s = 0;
% --------------------------------------------------------------
function u0 = icfun(x)
global beta_p100 beta_p0 c_smaxp I
if I < 0
    u0 = beta_p100*c_smaxp; % Discharge initial concentration (t=0)
end
if I > 0
    u0 = beta_p0*c_smaxp;   % Charge initial concentration (t=0)
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I t_corte
F = 96478;
a = 3*0.5/1e-4;
L = 36.4e-4;
D_s = 3.7e-12;

pr = 0;
if  (t <= t_corte)
    pr = (I)./(F*a*L);
end

pl = 0;
ql = 1;
qr = D_s;