function [c_sp] = concentration_p(r,t)
global c_smaxp beta_p100 beta_p0

m = 2;
x=r;

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
u = sol(:,:,1); % Extract the first solution component as u

c_sp = u;

% Calculation of c_s average
for i = 1:1:length(t)
    Int = 0;
    for j = 1:1:(length(r)-1)
        Int = Int + (1/2)*(r(j)^2*u(i,j) + r(j+1)^2*u(i,j+1))*(r(j+1)-r(j));
    end
    c_sp_prom(i) = Int.*(3/r(end)^3);
end

% Surface plot
figure(1)
surf(x,t,u) 
title('Concentration electrode^{+}')
xlabel('Distance, r [cm]')
ylabel('Time, t [s]')
zlabel('Concentration^+, c_{s} [mol*m^{-3}]')

% Average concentration
figure(2)
plot(t,c_sp_prom)
hold on
plot(t,c_smaxp*beta_p0*ones(1,length(t)),'r')
plot (t,c_smaxp*beta_p100*ones(1,length(t)))
plot(t,u(:,end),'g')
plot(t,u(:,1),'k')
title('Average concentration electrode^{+}')
xlabel('Time, t [s]')
ylabel('c_{sp,prom} [mol*s^{-1}]')
legend('c_{sp,prom}', 'max limit', 'min limit', 'c_{sp}(t,R_p)', 'c_{sp}(t,0)')

% Solution profile.
figure(5)
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
global beta_p100 c_smaxp

u0 = beta_p100*c_smaxp; % Initial concentration (t=0)
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I t_corte
F = 96478;
a = 3*0.5/1e-4;
L = 36.4e-4;
D_s = 3.7e-12;

pr = 0;
if  (t < t_corte)
    pr = (I)./(F*a*L);
end

pl = 0;
ql = 1;
qr = D_s;