function [c_sn] = concentration_n(r,t)
global c_smaxn beta_n100 beta_n0

m = 2;
x=r;

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
u = sol(:,:,1); % Extract the first solution component as u

c_sn = u;

% Calculation of c_s average
for i = 1:1:length(t)
    Int = 0;
    for j = 1:1:(length(r)-1)
        Int = Int + (1/2)*(r(j)^2*u(i,j) + r(j+1)^2*u(i,j+1))*(r(j+1)-r(j));
    end
    c_sn_prom(i) = Int.*(3/r(end)^3);
end

% Surface plot
figure(3)
surf(x,t,u) 
title('Concentration electrode^{-}')
xlabel('Distance, r [cm]')
ylabel('Time, t [s]')
zlabel('Concentration^-, c_{s} [mol*m^{-3}]')

% Average concentration
figure(4)
plot(t,c_sn_prom)
hold on
plot (t,c_smaxn*beta_n100*ones(1,length(t)),'r')
plot(t,c_smaxn*beta_n0*ones(1,length(t)))
plot(t,u(:,end),'g')
plot(t,u(:,1),'k')
title('Average concentration electrode^{-}')
xlabel('Time, t [s]')
ylabel('c_{sn,prom} [mol*s^{-1}]')
legend('c_{sn,prom}', 'max limit', 'min limit', 'c_{sn}(t,R_p)', 'c_{sn}(t,0)')

% Solution profile
figure(5)
subplot(2,1,2);
plot(t,u(:,end))
title('Concentration electrode^{-} at R_{p}')
xlabel('Time, t [s]')
ylabel('c_{s}(R_{p},t) [mol*cm^{-3}]')

figure(10)
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