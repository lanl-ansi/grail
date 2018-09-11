function [rho]=p_to_rho(p,b1,b2,psc)
rho=p.*(b1+b2*psc*p);