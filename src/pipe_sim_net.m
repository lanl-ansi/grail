function [tt,pp,qq]=pipe_sim_net(par)
%everything nondimensionalized for simulation

FN=par.FN; M=par.M; NE=par.NE;

%initial conditions
x0=[par.ppp0;par.qqq0];

%most of this is obsolete
part=par;
Tsc=par.Tsc; T=par.Ts*par.Tsc;
part.dout=@(t) par.dout(mod(t*Tsc,T));
part.cfun=@(t) par.cfun(mod(t*Tsc,T));
part.sfun=@(t) par.sfun(mod(t*Tsc,T)); 
part.dsfun=@(t) par.dsfun(mod(t*Tsc,T)); 

Tspan=par.Tgridnd;
par.JS=abs(sign(pipe_net_f_imp_jac(1,rand(FN+NE,1),rand(FN+NE,1),par)));
opts=odeset('RelTol',part.rtol,'AbsTol',part.atol,...
    'Jacobian',@(t,x,dx) pipe_net_f_imp_jac(t,x,dx,part),'Jpattern',par.JS,'Vectorized','on');
[tt,xx] = ode15i(@(t,x,dx) pipe_net_f_imp(t,x,dx,part),Tspan,x0,pipe_net_f(0,x0,part),opts);
if(tt(end)<Tspan(end)) 
    disp(['error 2: t_end=' num2str(tt(end))]), plot(xx), return;
end
pp=xx(:,1:FN); qq=xx(:,FN+1:M);
tt=tt*par.Tsc; 
pp=xx(:,1:FN); qq=xx(:,FN+1:M);

%explicit ODE dynamics function
function [dx]=pipe_net_f(t,x,par)
Dk=par.Dk; lamk=par.lamk; ML=par.ML; Ma=par.Ma;
Ad=par.Ad; Am=par.Am; Xs=par.Xs; Adp=par.Adp; Amn=par.Amn;
%dsfun=par.dsfun;
FN=par.FN; NE=par.NE; M=par.M;
cn=par.comp_pos(:,1); cl=par.comp_pos(:,2);
dnodes=par.dnodes; snodes=par.snodes;

p=x(1:FN); q=x(FN+1:M);     %the variables
s=par.sfun(t);      %this is slack node density
ds=par.dsfun(t);    %derivative of slack node density
comps=par.cfun(t);  %compression values at time t

if(par.doZ==1), b1=par.b1; b2=par.b2; psc=par.psc; 
s=p_to_rho(s,b1,b2,psc); end

Amj=Am; Amj(sub2ind(size(Amj),cn,cl))=-comps;
Adj=Amj(dnodes,:); Asj=Amj(snodes,:);
Amnj=Amn; Amnj(sub2ind(size(Amnj),cn,cl))=-comps;
Adnj=Amnj(dnodes,:);
if(par.doZ==1), AsAdp2=(-rho_to_p(-Asj'*s,b1,b2,psc)+rho_to_p(Adp'*p,b1,b2,psc)-rho_to_p(-Adnj'*p,b1,b2,psc)); 
else AsAdp2=(Asj'*s+Adj'*p); end
dx=(abs(Ad)*Xs*ML*abs(Adj'))\(Ad*Xs*q-par.dout(t)-abs(Ad)*Xs*ML*abs(Asj')*ds);                                  %density dynamics
dx=[dx;Ma*AsAdp2-q.*abs(q)./(abs(Asj')*s+abs(Adj')*p).*lamk./Dk*par.R];   %flux dynamics

%implicit DAE dynamics function
function [ff]=pipe_net_f_imp(t,x,dx,par)
Dk=par.Dk; lamk=par.lamk; ML=par.ML; Ma=par.Ma;
Xs=diag(sparse(Dk.^2/4*pi));
Ad=par.Ad; Am=par.Am; Adp=par.Adp; Amn=par.Amn;
FN=par.FN; NE=par.NE; M=par.M;
cn=par.comp_pos(:,1); cl=par.comp_pos(:,2);
dnodes=par.dnodes; snodes=par.snodes;

p=x(1:FN); q=x(FN+1:M);     %the variables
s=par.sfun(t);      %this is slack bus density
ds=par.dsfun(t);    %derivative of slack node density
dp=dx(1:FN); dq=dx(FN+1:M);     %the derivatives
comps=par.cfun(t);  %compression values at time t

if(par.doZ==1), b1=par.b1; b2=par.b2; psc=par.psc; 
s=p_to_rho(s,b1,b2,psc); end

Amj=Am; Amj(sub2ind(size(Amj),cn,cl))=-comps;
Adj=Amj(dnodes,:); Asj=Amj(snodes,:);
Amnj=Amn; Amnj(sub2ind(size(Amnj),cn,cl))=-comps;
Adnj=Amnj(dnodes,:);
if(par.doZ==1), AsAdp2=(-rho_to_p(-Asj'*s,b1,b2,psc)+rho_to_p(Adp'*p,b1,b2,psc)-rho_to_p(-Adnj'*p,b1,b2,psc)); 
else AsAdp2=(Asj'*s+Adj'*p); end
ff=(abs(Ad)*Xs*ML*abs(Adj'))*dp-(Ad*Xs*q-par.dout(t)-abs(Ad)*Xs*ML*abs(Asj')*ds);   %density dynamics
ff=[ff;(abs(Asj')*s+abs(Adj')*p).*dq-Ma*(AsAdp2.*(abs(Asj')*s+abs(Adj')*p))+q.*abs(q).*lamk./Dk*par.R];   %flux dynamics

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

if(par.doZ==1), b1=par.b1; b2=par.b2; psc=par.psc; 
s=p_to_rho(s,b1,b2,psc); end

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