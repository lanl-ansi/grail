function [cx]=pipe_grad_base_s(x,par)

%index parameters
N1=par.N1; N=par.N; M=par.M; C=par.C; FN=par.FN; NE=par.NE; GN=par.GN;
Ts=par.Ts; w=par.w; D=-par.D*2/Ts; prd=par.prd; prslack=par.prslack;
clinks=par.comp_pos(:,2); slinks=par.comp_pos(par.spos,2);
%problem parameters
m=par.mpow; xsc=par.xs(clinks); ew=par.econweight;
%variables
if(par.use_init_state==0)
    q=reshape(x(N1*FN+1:N1*M),NE,N1);
    comps=reshape(x(N1*M+1:N1*M+N1*C),C,N1);
    fd=reshape(x(N1*(M+C)+1:N1*(M+C+GN)),GN,N1);
elseif(par.use_init_state==1)
    st_ind=FN;
    %pp_init=par.state_init(st_ind+1:st_ind+FN); st_ind=st_ind+FN;
    qq_init=par.state_init(st_ind+1:st_ind+NE); st_ind=st_ind+NE;
    cc_init=par.state_init(st_ind+1:st_ind+C); st_ind=st_ind+C;
    dd_init=par.state_init(st_ind+1:st_ind+GN); st_ind=st_ind+GN;
    %p=[pp_init reshape(x(1:N*FN),FN,N)]; 
    q=[qq_init reshape(x(N*FN+1:N*M),NE,N)];
    comps=[cc_init reshape(x(N*M+1:N*M+N*C),C,N)];
    fd=[dd_init reshape(x(N*(M+C)+1:N*(M+C+GN)),GN,N)];
end
qcomp=q(clinks,:); fs=q(slinks,:);

% %restrict objective to time interval of interest only
% if(par.extension>0)
%     Nf=N1-par.extension;
%     D=D(1:Nf,1:Nf);
%     qcomp=qcomp(:,1:Nf); comps=comps(:,1:Nf); w=w(1:Nf);
%     fd=fd(:,1:Nf); prd=prd(:,1:Nf); prslack=prslack(:,1:Nf); fs=fs(:,1:Nf);
% end

%initialize, 
%cxeff=sparse(zeros(N1*(M+C+GN),1));
%cxecon=sparse(zeros(N1*(M+C+GN),1));
cxeff=zeros(N1*(M+C+GN),1);
cxecon=zeros(N1*(M+C+GN),1);
for j=1:C
    %dcdq
    cxeff(N1*FN+clinks(j):NE:N1*FN+NE*(N1-1)+clinks(j))=xsc(j)*sign(qcomp(j,:)').*((comps(j,:)').^(m)-1).*w;
    %dcdc
    cxeff(N1*M+j:C:N1*M+C*(N1-1)+j)=m*abs(qcomp(j,:)'*xsc(j)).*((comps(j,:)').^(m-1)).*w;
end
cxeff=cxeff; cxsmc=[zeros(N1*M,1);par.smsc*reshape((2*(D'*diag(w)*D)*comps')',C*N1,1);zeros(N1*GN,1)];
for j=1:length(slinks)
    %dedfs
    cxecon(N1*FN+slinks(j):NE:N1*FN+NE*(N1-1)+slinks(j))=(xsc(j)*prslack(j,:)').*w;
end
for j=1:GN
    %dedfd
    cxecon(N1*(M+C)+j:GN:N1*(M+C)+GN*(N1-1)+j)=-prd(j,:)'.*w;
end
cxecon=cxecon+[zeros(N1*(M+C),1);par.smsd*reshape((2*(D'*diag(w)*D)*fd')',GN*N1,1)];
cx=(Ts/2)*(cxeff*(1-ew)+cxecon*ew+cxsmc)/par.objsc;

if(par.use_init_state==1)
    col_elim=[[1:FN]'; [N1*FN+1:N1*FN+NE]'; [N1*M+1:N1*M+C]'; [N1*(M+C)+1:N1*(M+C)+GN]'];
    cx(col_elim)=[];
end