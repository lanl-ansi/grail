function [dp]=rho_to_p_diff(rho,b1,b2,psc)
dp=rho./sqrt(b1^2+4*b2*psc*rho);