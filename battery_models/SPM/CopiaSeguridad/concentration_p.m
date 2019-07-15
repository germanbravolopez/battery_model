function [c_sp] = concentration_p(r,t,I);

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
plot(x,u(end,:))
title('Concentration electrode^{+} at t_{final}')
xlabel('Distance, r [cm]')
ylabel('c_{s}(r,t_{final}) [mol*cm^{-3}]')
% --------------------------------------------------------------
function [c,f,s] = pdefun(x,t,u,DuDx)
D_s = 3.7e-12;

c = 1;
f = D_s*DuDx;
s = 0;
% --------------------------------------------------------------
function u0 = icfun(x)
global beta_p100 c_smaxp

u0 = beta_p100*c_smaxp; % Concentración inicial (t=0)
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I c_smaxp beta_p100 beta_p0 %cont cont1
F = 96478;
a = 3*0.5/1e-4;
L = 36.4e-4;
D_s = 3.7e-12;

if  (ur < c_smaxp*beta_p100) || (ur > c_smaxp*beta_p0) %|| (cont1 >= cont)
    I = 0;
end
% cont1 = cont1 + 1;

pl = 0;
ql = 1/D_s;
pr = (I)./(F*a*L);
qr = 1;