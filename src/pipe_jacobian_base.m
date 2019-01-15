function [J]=pipe_jacobian_base(x,par)

%index parameters
N1=par.N1; N=par.N; M=par.M; C=par.C; FN=par.FN; NE=par.NE; GN=par.GN;
Ts=par.Ts; w=par.w; D=-par.D*2/Ts; dnodes=par.dnodes; snodes=par.snodes;
clinks=par.comp_pos(:,2); slinks=par.comp_pos(par.spos,2);
%problem parameters
m=par.mpow; Xs=par.Xs; xs=par.xs;
Dk=par.Dk; lamk=par.lamk; ML=par.ML; Ma=par.Ma;
Am=par.Am; Ad=par.Ad; As=par.As; d=par.d; R=par.R; Adp=par.Adp; Adn=par.Adn; Amn=par.Amn;
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
cpos=par.comp_pos; spos=par.spos; dpos=par.dpos; ppos=par.ppos; sppos=par.sppos;
s=par.pslack;
%scf=find(ismember(par.snodes,cpos(spos)));
%s(scf,:)=comps(spos,:).*s(scf,:);
d(par.guniqueind,:)=d(par.guniqueind,:)+par.gtod*fd;
qcomp=q(clinks,:);

if(par.doZ==1), b1=par.b1; b2=par.b2; psc=par.psc; end

J=sparse(N1*M,N1*(M+C+GN));        %allocate space
%define indices
p1=1;p2=N1*FN;q1=p2+1;q2=p2+N1*NE;c1=q2+1;c2=q2+C*N1;
%dynamic constraints
for j=1:N1
    Amj=Am; 
    Amj(sub2ind(size(Amj),cpos(:,1),cpos(:,2)))=-comps(:,j);
    Adj=Amj(dnodes,:); Asj=Amj(snodes,:); 
    Amnj=Amn; Amnj(sub2ind(size(Amnj),cpos(:,1),cpos(:,2)))=-comps(:,j);
    Adnj=Amnj(par.dnodes,:);
    pDjt=p*(D(j,:))'; qDjt=q*(D(j,:))'; sDjt=s*(D(j,:))'; 
    AsAdpj=(abs(Asj')*s(:,j)+abs(Adj')*p(:,j));
    if(par.doZ==1) 
        AsAdp2j=(-rho_to_p(-Asj'*s(:,j),b1,b2,psc)+rho_to_p(Adp'*p(:,j),b1,b2,psc)-rho_to_p(-Adnj'*p(:,j),b1,b2,psc));
        dAdjdp2j=diag(sparse(rho_to_p_diff(Adp'*p(:,j),b1,b2,psc)))*Adp'...
            +diag(sparse(rho_to_p_diff(-Adnj'*p(:,j),b1,b2,psc)))*Adnj';
    else
        AsAdp2j=(Asj'*s(:,j)+Adj'*p(:,j));
        dAdjdp2j=Adj';
    end
    %ddp/dp
    J(FN*(j-1)+1:FN*j,p1:p2)=kron(D(j,:),-(abs(Ad)*Xs*ML*abs(Ad')));
    %ddp/dq
    J(FN*(j-1)+1:FN*j,p2+NE*(j-1)+1:p2+NE*j)=Ad*Xs;
    %ddp/dfd
    %J(FN*(j-1)+1:FN*j,c2+GN*(j-1)+1:c2+GN*j)=-sparse(par.gallind,[1:GN]',ones(GN,1),FN,GN);
    J(FN*(j-1)+1:FN*j,c2+GN*(j-1)+1:c2+GN*j)=-sparse(par.gallind,[1:GN]',par.gallsign,FN,GN);
    %ddq/dp
    J(p2+NE*(j-1)+1:p2+NE*j,FN*(j-1)+1:FN*j)=Ma*(diag(sparse(AsAdpj))*dAdjdp2j+...
        diag(sparse(AsAdp2j))*abs(Adj'))-diag(sparse(qDjt))*abs(Adj');
    %ddq/dq
    J(p2+NE*(j-1)+1:p2+NE*j,q1:q2)=kron(D(j,:),-diag(sparse(AsAdpj)));
    J(p2+NE*(j-1)+1:p2+NE*j,p2+NE*(j-1)+1:p2+NE*j)=J(p2+NE*(j-1)+1:p2+NE*j,p2+NE*(j-1)+1:p2+NE*j)...
        +R*gxfun(q(:,j),lamk,Dk);
    %----------------
    ddPos=sparse(cpos(dpos,2),[1:length(dpos)],pDjt(ppos),NE,length(dpos));
    dPos=sparse(cpos(dpos,2),[1:length(dpos)],-p(ppos,j),NE,length(dpos));
    if(par.doZ==1) 
        dPosz=rho_to_p_diff(comps(dpos,j).*p(ppos,j),b1,b2,psc).*p(ppos,j);
        dPosz1=sparse(cpos(dpos,2),[1:length(dpos)],-dPosz,NE,length(dpos));
    else
        dPosz1=dPos;
    end    
    %ddp/da
    J(FN*(j-1)+1:FN*j,q2+C*(j-1)+dpos)=-abs(Ad)*Xs*ML*ddPos;
    %ddq/da
    J(p2+NE*(j-1)+1:p2+NE*j,q2+C*(j-1)+dpos)= Ma*(diag(sparse((AsAdpj)))*dPosz1-diag(sparse((AsAdp2j)))*dPos)...
        +diag(sparse(qDjt))*dPos;
    %----------------
    ssPos=sparse(cpos(spos,2),[1:length(spos)],sDjt,NE,length(spos));
    sPos=sparse(cpos(spos,2),[1:length(spos)],-s(:,j),NE,length(spos));
    if(par.doZ==1) 
        sPosz=rho_to_p_diff(comps(spos,j).*s(:,j),b1,b2,psc).*s(:,j);
        sPosz1=sparse(cpos(spos,2),[1:length(spos)],-sPosz,NE,length(spos));
    else
        sPosz1=sPos;
    end
    %sPosz=sparse(cpos(spos,2),[1:length(spos)],-dpsz(:,j).*s(:,j)./comps(spos,j),NE,length(spos));
    %ddp/ds
    J(FN*(j-1)+1:FN*j,q2+C*(j-1)+spos)= -abs(Ad)*Xs*ML*ssPos;
    %ddq/ds
    J(p2+NE*(j-1)+1:p2+NE*j,q2+C*(j-1)+spos)= Ma*(diag(sparse((AsAdpj)))*sPosz1-diag(sparse((AsAdp2j)))*sPos)...
        +diag(sparse(qDjt))*sPos;
end

if(par.use_init_state==0)
    %discharge pressure (density) constraints
    Jineq=sparse(C*N1,N1*(M+C+GN));
    for j=1:length(spos)
        II=[spos(j):C:(N1-1)*C+spos(j)]';
        JJ=[q2+spos(j):C:q2+(N1-1)*C+spos(j)]';
        %SS=par.p_min_nd(cpos(spos(j),1))*ones(N1,1);
        SS=par.pslack(sppos(j),:)';
        Jineq=Jineq+sparse(II,JJ,SS,C*N1,N1*(M+C+GN));
    end
    for j=1:length(dpos)
        II=[[dpos(j):C:(N1-1)*C+dpos(j)]';[dpos(j):C:(N1-1)*C+dpos(j)]'];
        JJ=[[q2+dpos(j):C:q2+(N1-1)*C+dpos(j)]';[ppos(j):FN:(N1-1)*FN+ppos(j)]'];
        SS=[p(ppos(j),:)';comps(dpos(j),:)'];
        Jineq=Jineq+sparse(II,JJ,SS,C*N1,N1*(M+C+GN));
    end
    J=[J;Jineq];

    %max power constraint
    Jineq2=sparse(C*N1,N1*(M+C+GN));
    for j=1:C
        II=[[j:C:(N1-1)*C+j]'; [j:C:(N1-1)*C+j]'];
        JJ=[[N1*FN+clinks(j):NE:N1*FN+NE*(N1-1)+clinks(j)]';[N1*M+j:C:N1*M+C*(N1-1)+j]'];
        SS=[[xs(clinks(j))*sign(qcomp(j,:)').*((comps(j,:)').^(m)-1)];...
        [m*abs(qcomp(j,:)'*xs(clinks(j))).*((comps(j,:)').^(m-1))]];
    Jineq2=Jineq2+sparse(II,JJ,SS,N1*C,N1*(M+C+GN));
    end
    J=[J;Jineq2*par.hplsc];
elseif(par.use_init_state==1)
        %discharge pressure (density) constraints
    Jineq=sparse(C*N,N1*(M+C+GN));
    for j=1:length(spos)
        II=[spos(j):C:(N-1)*C+spos(j)]';
        JJ=[q2+spos(j)+C:C:q2+(N1-1)*C+spos(j)]';
        %SS=par.p_min_nd(cpos(spos(j),1))*ones(N1,1);
        SS=par.pslack(sppos(j),2:N1)';
        Jineq=Jineq+sparse(II,JJ,SS,C*N,N1*(M+C+GN));
    end
    for j=1:length(dpos)
        II=[[dpos(j):C:(N-1)*C+dpos(j)]';[dpos(j):C:(N-1)*C+dpos(j)]'];
        JJ=[[q2+dpos(j)+C:C:q2+(N1-1)*C+dpos(j)]';[ppos(j)+FN:FN:(N1-1)*FN+ppos(j)]'];
        SS=[p(ppos(j),2:N1)';comps(dpos(j),2:N1)'];
        Jineq=Jineq+sparse(II,JJ,SS,C*N,N1*(M+C+GN));
    end
    J=[J;Jineq];

    %max power constraint
    Jineq2=sparse(C*N,N1*(M+C+GN));
    for j=1:C
        II=[[j:C:(N-1)*C+j]'; [j:C:(N-1)*C+j]'];
        JJ=[[N1*FN+clinks(j)+NE:NE:N1*FN+NE*(N1-1)+clinks(j)]';[N1*M+j+C:C:N1*M+C*(N1-1)+j]'];
        SS=[[xs(clinks(j))*sign(qcomp(j,2:N1)').*((comps(j,2:N1)').^(m)-1)];...
        [m*abs(qcomp(j,2:N1)'*xs(clinks(j))).*((comps(j,2:N1)').^(m-1))]];
    Jineq2=Jineq2+sparse(II,JJ,SS,N*C,N1*(M+C+GN));
    end
    J=[J;Jineq2*par.hplsc];
end

if(par.int_flow_const==1)
    J=[J;par.Aineq];
end

%if initial state constraint, eliminate fixed variable columns
if(par.use_init_state==1)
    row_elim=[[1:FN]'; [N1*FN+1:N1*FN+NE]'];
    col_elim=[[1:FN]'; [N1*FN+1:N1*FN+NE]'; [N1*M+1:N1*M+C]'; [N1*(M+C)+1:N1*(M+C)+GN]'];
    J(:,col_elim)=[]; J(row_elim,:)=[];
end

%add constant discharge interval constraints
%J=[par.Aeq;J];

function [g]=gfun(x,lamk,Dk)
g=-(lamk./Dk).*x.*abs(x);

function [gx]=gxfun(x,lamk,Dk)
gx=diag(sparse(-2*(lamk./Dk).*abs(x)));

function [p]=rho_to_p(rho,b1,b2,psc)
p=(-b1+sqrt(b1^2+4*(b2*psc)*rho))/(2*b2*psc);

function [dp]=rho_to_p_diff(rho,b1,b2,psc)
dp=1./sqrt(b1^2+4*(b2*psc)*rho);

