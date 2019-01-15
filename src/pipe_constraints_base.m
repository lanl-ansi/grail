function [k]=pipe_constraints_base(x,par)

%index parameters
N1=par.N1; N=par.N; M=par.M; C=par.C; FN=par.FN; NE=par.NE; GN=par.GN;
Ts=par.Ts; w=par.w; D=-par.D*2/Ts; 
clinks=par.comp_pos(:,2); slinks=par.comp_pos(par.spos,2);
%problem parameters
m=par.mpow; Xs=par.Xs; 
Dk=par.Dk; lamk=par.lamk; ML=par.ML; Ma=par.Ma; Adp=par.Adp; Adn=par.Adn; Amn=par.Amn;
Am=par.Am; Ad=par.Ad; As=par.As; d=par.d; R=par.R;
%variables
if(par.use_init_state==0)
    p=reshape(x(1:N1*FN),FN,N1); q=reshape(x(N1*FN+1:N1*M),NE,N1);
    comps=reshape(x(N1*M+1:N1*M+N1*C),C,N1);
    fd=reshape(x(N1*(M+C)+1:N1*(M+C+GN)),GN,N1);
elseif(par.use_init_state==1)
    st_ind=0;
    pp_init=par.state_init(st_ind+1:st_ind+FN); st_ind=st_ind+FN;
    qq_init=par.state_init(st_ind+1:st_ind+NE); st_ind=st_ind+NE;
    cc_init=par.state_init(st_ind+1:st_ind+C); st_ind=st_ind+C;
    dd_init=par.state_init(st_ind+1:st_ind+GN); st_ind=st_ind+GN;
    p=[pp_init reshape(x(1:N*FN),FN,N)]; 
    q=[qq_init reshape(x(N*FN+1:N*M),NE,N)];
    comps=[cc_init reshape(x(N*M+1:N*M+N*C),C,N)];
    fd=[dd_init reshape(x(N*(M+C)+1:N*(M+C+GN)),GN,N)];
end
cpos=par.comp_pos; spos=par.spos; dpos=par.dpos; ppos=par.ppos;
s=par.pslack;
%scf=find(ismember(par.snodes,cpos(spos)));
%s(scf,:)=comps(spos,:).*s(scf,:);
d(par.guniqueind,:)=d(par.guniqueind,:)+par.gtod*fd;
qcomp=q(clinks,:);

if(par.doZ==1), b1=par.b1; b2=par.b2; psc=par.psc; end

if(par.use_init_state==0), k=zeros(N1*(M+2*C),1); %k=sparse(N1*(M+2*C),1);
elseif(par.use_init_state==1), k=zeros(N1*M+2*N*C,1); end %k=sparse(N1*M+2*N*C,1); end
for j=1:N1
   Amj=Am; Amj(sub2ind(size(Amj),cpos(:,1),cpos(:,2)))=-comps(:,j);
   Adj=Amj(par.dnodes,:); Asj=Amj(par.snodes,:);
   Amnj=Amn; Amnj(sub2ind(size(Amnj),cpos(:,1),cpos(:,2)))=-comps(:,j);
   Adnj=Amnj(par.dnodes,:);
   if(par.doZ==1), AsAdp2j=(-rho_to_p(-Asj'*s(:,j),b1,b2,psc)+rho_to_p(Adp'*p(:,j),b1,b2,psc)-rho_to_p(-Adnj'*p(:,j),b1,b2,psc)); 
   else AsAdp2j=(Asj'*s(:,j)+Adj'*p(:,j)); end
  k(FN*(j-1)+1:FN*j)=...
       (Ad*Xs*q(:,j)-d(:,j)-abs(Ad)*Xs*ML*abs(Asj')*s*(D(j,:))')-(abs(Ad)*Xs*ML*abs(Adj'))*p*(D(j,:))';
  k(N1*FN+NE*(j-1)+1:N1*FN+NE*j)=Ma*AsAdp2j.*(abs(Asj')*s(:,j)+abs(Adj')*p(:,j))+...
       gfun(q(:,j),lamk,Dk)*R-(q*(D(j,:))').*(abs(Asj')*s(:,j)+abs(Adj')*p(:,j));
end

if(par.use_init_state==0)
    %discharge pressure constraint
    rhosuction=sparse(C,N1);
    rhosuction(spos,:)=par.pslack;
    rhosuction(dpos,:)=p(ppos,:);
    if(par.doZ==1)
    k(N1*M+1:N1*(M+C),1)=-kron(ones(N1,1),p_to_rho(par.p_max_nd(par.comp_pos(:,1)),b1,b2,psc))...
        +reshape((rhosuction.*comps),C*N1,1);
    else
        k(N1*M+1:N1*(M+C),1)=-kron(ones(N1,1),par.p_max_nd(par.comp_pos(:,1)))...
        +reshape((rhosuction.*comps),C*N1,1);
    end
    %max power constraint
    k(N1*(M+C)+1:N1*(M+2*C),1)=(-kron(ones(N1,1),par.boost_pow_max_nd)+...
        reshape((Xs(cpos(:,2),cpos(:,2))*abs(qcomp)).*((comps).^(m)-1),C*N1,1))*par.hplsc;
elseif(par.use_init_state==1)
    rhosuction=sparse(C,N);
    rhosuction(spos,:)=par.pslack(:,2:N1);
    rhosuction(dpos,:)=p(ppos,2:N1);
    if(par.doZ==1)
    k(N1*M+1:N1*M+N*C,1)=-kron(ones(N,1),p_to_rho(par.p_max_nd(par.comp_pos(:,1)),b1,b2,psc))...
        +reshape((rhosuction.*comps(:,2:N1)),C*N,1);
    else
        k(N1*M+1:N1*M+N*C,1)=-kron(ones(N,1),par.p_max_nd(par.comp_pos(:,1)))...
        +reshape((rhosuction.*comps(:,2:N1)),C*N,1);
    end
    %max power constraint
    k(N1*M+N*C+1:N1*M+2*N*C,1)=(-kron(ones(N,1),par.boost_pow_max_nd)+...
        reshape((Xs(cpos(:,2),cpos(:,2))*abs(qcomp(:,2:N1))).*((comps(:,2:N1)).^(m)-1),C*N,1))*par.hplsc;
end

if(par.int_flow_const==1)
    k=[k;par.Aineq*x];
end

if(par.use_init_state==1)
    row_elim=[[1:FN]'; [N1*FN+1:N1*FN+NE]'];
    k(row_elim)=[];
end

%add constant discharge interval constraints
%k=[par.Aeq*x-par.Beq;k];

%par.cintl=par.N1/par.compints;

function [g]=gfun(x,lamk,Dk)
g=-(lamk./Dk).*x.*abs(x);

function [p]=rho_to_p(rho,b1,b2,psc)
p=(-b1+sqrt(b1^2+4*b2*psc*rho))/(2*b2*psc);

function [rho]=p_to_rho(p,b1,b2,psc)
rho=p.*(b1+b2*psc*p);