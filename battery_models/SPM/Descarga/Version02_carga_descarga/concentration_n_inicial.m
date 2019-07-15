function [u] = concentration_n_inicial(r,t)

m = 2;
x = r;

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);

figure(2)
surf(x,t,u)

% figure(4)
% subplot(2,1,2);
% plot(t,c_sn_prom)
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
global beta_n100 beta_n0 c_smaxn I
if I < 0
    u0 = beta_n100*c_smaxn; % Discharge initial concentration (t=0)
end
if I > 0
    u0 = beta_n0*c_smaxn;   % Charge initial concentration (t=0)
end
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I beta_n0 beta_n100 c_smaxn
F = 96478;
a = 3*0.58/1e-4;
L = 50e-4;
D_s = 2e-12;

if  (ur < c_smaxn*beta_n0) || (ur > c_smaxn*beta_n100) 
    I = 0;
end

pl = 0;
ql = 1;
pr = -(I)./(F*a*L);
qr = D_s;
end