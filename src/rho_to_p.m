function [p]=rho_to_p(rho,b1,b2)
p=(-b1+sqrt(b1^2+4*b2*rho))/(2*b2);