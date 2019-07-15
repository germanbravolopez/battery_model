function [c_sn_prom] = concentration_n_inicial(r,t)
global c_smaxn beta_n100 beta_n0 

m = 2;
x = r; % Copy the input radius 

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
u = sol(:,:,1); % Extract the first solution component as u

% Calculation of c_s average
for i = 1:1:length(t)
    Int = 0;
    for j = 1:1:(length(r)-1)
        Int = Int + (1/2)*(r(j)^2*u(i,j) + r(j+1)^2*u(i,j+1))*(r(j+1)-r(j));
    end
    c_sn_prom(i) = Int.*(3/r(end)^3);
    if (c_sn_prom(i) < c_smaxn*beta_n0) || (c_sn_prom(i) > c_smaxn*beta_n100)
        c_sn_prom(i) = 0.02;
    end
end

% figure(3)
% surf(x,t,u) 
% title('Concentration electrode^{-}')
% xlabel('Distance, r [cm]')
% ylabel('Time, t [s]')
% zlabel('Concentration^-, c_{s} [mol*m^{-3}]')
% 
% figure(4)
% plot(t,c_sn_prom)
% hold on
% plot (t,c_smaxn*beta_n100*ones(1,length(t)),'r')
% plot(t,c_smaxn*beta_n0*ones(1,length(t)),'r')
% plot(t,u(:,end),'g')
% plot(t,u(:,1),'k')
% title('Average concentration electrode^{-}')
% xlabel('Time, t [s]')
% ylabel('c_{sn,prom} [mol*s^{-1}]')

end
% --------------------------------------------------------------
function [c,f,s] = pdefun(x,t,u,DuDx)
D_s = 2e-12;

c = 1/D_s;
f = DuDx;
s = 0;
end
% --------------------------------------------------------------
function u0 = icfun(x)
global beta_n100 c_smaxn

u0 = beta_n100*c_smaxn; % Initial concentration (t=0)
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I c_smaxn
F = 96478;
a = 3*0.58/1e-4;
L = 50e-4;
D_s = 2e-12;

if  (ur/c_smaxn > 1) || (ur/c_smaxn < 0) || (ul/c_smaxn > 1) || (ul/c_smaxn < 0) 
    I = 0;
end

pl = 0;
ql = 1;
pr = -(I)./(F*a*L);
qr = D_s;
end