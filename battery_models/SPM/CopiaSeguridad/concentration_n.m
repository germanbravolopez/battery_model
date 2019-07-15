function [c_sn] = concentration_n(r,t);

m = 2;
x=r;

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);
c_sn = u;

% A surface plot is often a good way to study a solution.
figure(2)
surf(x,t,u) 
title('Concentration electrode^{-}.')
xlabel('Distance, r [cm]')
ylabel('Time, t [s]')
zlabel('Concentration^-, c_{s} [mol*m^{-3}]')

% A solution profile can also be illuminating.
figure(3)
subplot(2,1,2);
plot(x,u(end,:))
title('Concentration electrode^{-} at t_{final}')
xlabel('Distance, r [cm]')
ylabel('c_{s}(r,t_{final}) [mol*cm^{-3}]')
% --------------------------------------------------------------
function [c,f,s] = pdefun(x,t,u,DuDx)
D_s = 2e-12;

c = 1;
f = D_s*DuDx;
s = 0;
% --------------------------------------------------------------
function u0 = icfun(x)
global beta_n100 c_smaxn

u0 = beta_n100*c_smaxn; % Concentración inicial (t=0)
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I c_smaxn beta_n100 beta_n0 %cont cont2
F = 96478;
a = 3*0.58/1e-4;
L = 50e-4;
D_s = 2e-12;

if  (ur < c_smaxn*beta_n0) || (ur > c_smaxn*beta_n100) %|| (cont2 >= cont)
    I = 0;
end
%cont2 = cont2 + 1;

pl = 0;
ql = 1/D_s;
pr = -(I)./(F*a*L);
qr = 1;