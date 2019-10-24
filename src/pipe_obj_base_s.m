function [c]=pipe_obj_base_s(x,par)

%index parameters
N1=par.N1; N=par.N; M=par.M; C=par.C; FN=par.FN; NE=par.NE; GN=par.GN;
Ts=par.Ts; w=par.w; D=-par.D*2/Ts; prd=par.prd; prslack=par.prslack;
clinks=par.comp_pos(:,2); slinks=par.comp_pos(par.spos,2);
%problem parameters
m=par.mpow; xs=par.xs; ew=par.econweight;
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

%efficiency objective
ceff=sum((diag(xs(clinks))*abs(qcomp)).*((comps).^(m)-1),1)*w;
csm=par.smsc*trace((D*comps')'*diag(sparse(w))*(D*comps')); %derivative penalty on comps
%economic objective
cecon=-sum(fd.*prd,1)*w+sum((diag(xs(slinks))*prslack).*fs,1)*w;
dsm=par.smsd*trace((D*fd')'*diag(sparse(w))*(D*fd'));   %derivative penalty on demands

%output objective
c=(Ts/2)*((ceff)*(1-ew)+(cecon+dsm)*ew+csm)/par.objsc;