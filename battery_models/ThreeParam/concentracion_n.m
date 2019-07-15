function c_s = concentracion_n(j_nn,R_p,c_sn_ini)
syms c_s(t) j_nn(t) R_p c_sn_ini
eqn = diff(c_s,t) == -3/R_p*j_nn;
cond = c_s(0) == c_sn_ini;
c_s(t) = dsolve(eqn,cond);
end
