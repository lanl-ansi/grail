function [rho]=p_to_rho_nd(p,b1,b2,psc)
rho=p.*(b1+b2*psc*p);
