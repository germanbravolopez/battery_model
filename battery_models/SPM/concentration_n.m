function [c_sn] = concentration_n(r,t)

m = 2;
x=r;

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
u = sol(:,:,1); % Extract the first solution component as u

c_sn = u;

% Surface plot
figure(2)
surf(x,t,u) 
title('Concentration electrode^{-}')
xlabel('Distance, r [cm]')
ylabel('Time, t [s]')
zlabel('Concentration^-, c_{s} [mol*m^{-3}]')

% Solution profile
figure(3)
subplot(2,1,2);
plot(t,u(:,end))
title('Concentration electrode^{-} at R_{p}')
xlabel('Time, t [s]')
ylabel('c_{s}(R_{p},t) [mol*cm^{-3}]')
% --------------------------------------------------------------
function [c,f,s] = pdefun(x,t,u,DuDx)
D_s = 2e-12;

c = 1/D_s;
f = DuDx;
s = 0;
% --------------------------------------------------------------
function u0 = icfun(x)
global beta_n100 c_smaxn

u0 = beta_n100*c_smaxn; % Initial concentration (t=0)
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I t_corte 
F = 96478;
a = 3*0.58/1e-4;
L = 50e-4;
D_s = 2e-12;

pr = 0;
if  (t < t_corte)
    pr = -(I)./(F*a*L);
end

pl = 0;
ql = 1;
qr = D_s;