function [p]=rho_to_p_nd(rho,b1,b2,psc)
p=(-b1+sqrt(b1^2+4*b2*psc*rho))/(2*b2*psc);