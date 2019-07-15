function [c_sp_prom] = concentration_p_inicial(r,t)
global c_smaxp beta_p100 beta_p0 I

m = 2;
x = r;

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);

% Average concentration c_s calculus
if I < 0
    for i = 1:1:length(t)
        Int = 0;
        for j = 1:1:(length(r)-1)
           Int = Int + (1/2)*(r(j)^2*u(i,j) + r(j+1)^2*u(i,j+1))*(r(j+1)-r(j));
            c_sp_prom(i) = Int.*(3/r(end)^3);
            if (c_sp_prom(i) < c_smaxp*beta_p100) || (c_sp_prom(i) > c_smaxp*beta_p0)
                c_sp_prom(i) = 0;
            end
        end
    end
end 
if I > 0
    for i = 1:1:length(t)
        Int = 0;
        for j = 1:1:(length(r)-1)
           Int = Int + (1/2)*(r(j)^2*u(i,j) + r(j+1)^2*u(i,j+1))*(r(j+1)-r(j));
            c_sp_prom(i) = Int.*(3/r(end)^3);
            if (c_sp_prom(i) < c_smaxp*beta_p100) || (c_sp_prom(i) > c_smaxp*beta_p0)
                c_sp_prom(i) = 0.03;
            end
        end
    end
end

figure(1)
surf(x,t,u)

figure(4)
subplot(2,1,1);
plot(t,c_sp_prom)
title('Average concentration electrode^{+}')
xlabel('Time, t [s]')
ylabel('c_{sp,prom} [mol*s^{-1}]')

end
% --------------------------------------------------------------
function [c,f,s] = pdefun(x,t,u,DuDx)
D_s = 3.7e-12;

c = 1/D_s;
f = DuDx;
s = 0;
end
% --------------------------------------------------------------
function u0 = icfun(x)
global beta_p100 beta_p0 c_smaxp I
if I < 0
    u0 = beta_p100*c_smaxp; % Discharge initial concentration (t=0)
end
if I > 0
    u0 = beta_p0*c_smaxp;   % Charge initial concentration (t=0)
end
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global I
F = 96478;
a = 3*0.5/1e-4;
L = 36.4e-4;
D_s = 3.7e-12;

if  (ur < 0) || (ur > 1) 
    I = 0;
end

pl = 0;
ql = 1;
pr = (I)./(F*a*L);
qr = D_s;
end