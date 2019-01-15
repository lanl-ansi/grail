%implicit DAE dynamics Jacobians
function [ffx,ffdx]=pipe_net_f_imp_jac(t,x,dx,par)
Dk=par.Dk; lamk=par.lamk; ML=par.ML; Ma=par.Ma;
Xs=par.Xs; Ad=par.Ad; Am=par.Am; Adp=par.Adp; Amn=par.Amn;
FN=par.FN; NE=par.NE; M=par.M;
cn=par.comp_pos(:,1); cl=par.comp_pos(:,2);
dnodes=par.dnodes; snodes=par.snodes;

p=x(1:FN); q=x(FN+1:M);     %the variables
s=par.sfun(t);     %this is slack bus density
ds=par.dsfun(t);    %derivative of slack node density
dp=dx(1:FN); dq=dx(FN+1:M);     %the derivatives
comps=par.cfun(t);  %compression values at time t

if(par.doZ==1), b1=par.b1; b2=par.b2; psc=par.psc; end

Amj=Am; Amj(sub2ind(size(Amj),cn,cl))=-comps;
Adj=Amj(dnodes,:); Asj=Amj(snodes,:);
Amnj=Amn; Amnj(sub2ind(size(Amnj),cn,cl))=-comps;
Adnj=Amnj(dnodes,:);
if(par.doZ==1), AsAdp2=(-rho_to_p(-Asj'*s,b1,b2,psc)+rho_to_p(Adp'*p,b1,b2,psc)-rho_to_p(-Adnj'*p,b1,b2,psc));
    dAdjdp2=diag(sparse(rho_to_p_diff(Adp'*p,b1,b2,psc)))*Adp'...
            +diag(sparse(rho_to_p_diff(-Adnj'*p,b1,b2,psc)))*Adnj';
else AsAdp2=(Asj'*s+Adj'*p); dAdjdp2=Adj'; end
ffx=[sparse(FN,FN) -Ad*Xs];     %[dfp/dp dfp/dq]
ffx=[ffx; [diag(sparse(dq))*abs(Adj')-Ma*(diag(sparse((abs(Asj')*s+abs(Adj')*p)))*dAdjdp2+diag(sparse(AsAdp2))*abs(Adj')) 2*diag(sparse(abs(q)).*lamk./Dk*par.R)]];   %[dfq/dp dfq/dq]
ffdx=[(abs(Ad)*Xs*ML*abs(Adj')) sparse(FN,NE); sparse(NE,FN) diag(sparse(abs(Asj')*s+abs(Adj')*p))];   %Jacobian w.r.t. dx

function [p]=rho_to_p(rho,b1,b2,psc)
p=(-b1+sqrt(b1^2+4*b2*psc*rho))/(2*b2*psc);

function [dp]=rho_to_p_diff(rho,b1,b2,psc)
dp=1./sqrt(b1^2+4*(b2*psc)*rho);

function [rho]=p_to_rho(p,b1,b2,psc)
rho=p.*(b1+b2*psc*p);