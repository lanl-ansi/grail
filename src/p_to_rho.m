function [rho]=p_to_rho(p,b1,b2,par)
rho=par.c.psc*p.*(b1+b2*par.c.psc*p)/(par.c.gasR*par.c.gasT);