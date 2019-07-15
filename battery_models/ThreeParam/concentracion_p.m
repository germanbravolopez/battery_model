function c_s = concentracion_p(j_np,R_p,c_sp_ini)
syms c_s(t) j_np(t) R_p c_sp_ini
eqn = diff(c_s,t) == -3/R_p*j_np;
cond = c_s(0) == c_sp_ini;
c_s(t) = dsolve(eqn,cond);
end
