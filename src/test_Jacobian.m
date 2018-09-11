[par.tr]=test_jac_setup(par.tr);

%% test gradient
k=.0000001;
rr=rand(size(par.tr.m.xu0));
cx=pipe_grad_base(rr,par.tr.m);
for j=1:length(rr)
    vt=zeros(size(rr)); vt(j)=k;
    cxt(j)=(pipe_obj_base(rr+vt,par.tr.m)-pipe_obj_base(rr,par.tr.m))/k;
end
norm(cx-cxt')
plot(abs([cx cxt']))
%plot(imag([cx cxt']))
%plot((cx-cxt')./(cxt+(cxt==0))')

%% test gradient
k=.000001;
rr=rand(size(par.tr.m.xu0));
cx=pipe_grad_base_s(rr,par.tr.m);
for j=1:length(rr)
    vt=zeros(size(rr)); vt(j)=k;
    cxt(j)=(pipe_obj_base_s(rr+vt,par.tr.m)-pipe_obj_base_s(rr,par.tr.m))/k;
end
norm(cx-cxt')
plot(abs([cx cxt']))
%plot(imag([cx cxt']))
%plot((cx-cxt')./(cxt+(cxt==0))')

%% test Jacobian
k=.000000001;
rr=rand(size(par.tr.m.xu0));
J1=pipe_jacobian_base(rr,par.tr.m);
J1t=sparse([]);
for j=1:length(rr)
    vt=zeros(size(rr)); vt(j)=k;
    Cip=pipe_constraints_base(rr+vt,par.tr.m); Ci=pipe_constraints_base(rr,par.tr.m);
    J1t(:,j)=sparse((Cip-Ci)/k);
    if(mod(j,100)==0) disp(num2str(j)), end
end

[norm(full(J1-J1t)) max(max(abs(J1-J1t)))]
%%
N1=par.tr.m.N1;M=par.tr.m.M;C=par.tr.m.C;FN=par.tr.m.FN;
notcomps=[1:N1*(M+C)];
powconst=[N1*(M+C)+1:N1*(M+C+FN)];
[norm(full(J1(comps,:)-J1t(comps,:))) max(max(abs(J1(comps,:)-J1t(comps,:))))]
[norm(full(J1(powconst,:)-J1t(powconst,:))) max(max(abs(J1(powconst,:)-J1t(powconst,:))))]

%%
%for j=1:length(xu0)
N1=par.tr.m.N1;M=par.tr.m.M;C=par.tr.m.C;FN=par.tr.m.FN;NE=par.tr.m.NE;PN=par.tr.m.PN; GN=par.tr.m.GN;
pindp=[1:N1*FN]; pindq=[N1*FN+1:N1*M];
%for j=1:(N1*M)
for j=N1*M+1:N1*(M+C)
%for j=N1*M+1:N1*(M+C)
%   plot([J1(j,:)-J1t(j,:)]), title(num2str(j)), pause
%   plot([J1(pindq,j)-J1t(pindq,j)]), title(num2str(j)), pause
    plot([J1(pindq,j)-J1t(pindq,j)]), title(num2str(j)), pause
%   plot((J1(pindq,j)-J1t(pindq,j))), title(num2str(j)), pause
%plot((J1(:,j)-J1t(:,j))./(J1t(:,j)+(J1t(:,j)==0))), title(num2str(j)), pause
end
%%
%for j=1:length(xu0)
N1=par.tr.m.N1;M=par.tr.m.M;C=par.tr.m.C;FN=par.tr.m.FN;
for j=(N1*(M+C)+1):(N1*(M+C+FN))
plot([J1(:,j) J1t(:,j)]), title(num2str(j)), %axis([1 length(xu0) -100 100]),
pause
end

%%
for j=1:length(par.tr.m.xu0)
plot(J1(:,j)-J1t(:,j)), title(num2str(j)), %axis([1 length(xu0) -100 100]),
pause
end

%%
for j=1:length(rr)
N1=par.tr.m.N1;M=par.tr.m.M;C=par.tr.m.C;FN=par.tr.m.FN;
%for j=(N1*(M)+1):(N1*(M+C))
plot((J1(j,:)-J1t(j,:))./(J1t(j,:)+(J1t(j,:)==0))), title(num2str(j)), pause
end
%%

[par.ss]=test_jac_setup(par.ss);


%% test gradient
k=.000000001;
rr=rand(size(par.ss.m.xu0));
cx=pipe_grad_base(rr,par.ss.m);
for j=1:length(rr)
    vt=zeros(size(rr)); vt(j)=k;
    cxt(j)=(pipe_obj_base(rr+vt,par.ss.m)-pipe_obj_base(rr,par.ss.m))/k;
end
norm(cx-cxt')
plot(abs([cx cxt']))
%plot(imag([cx cxt']))
%plot((cx-cxt')./(cxt+(cxt==0))')

%% test gradient
k=.000000001;
rr=rand(size(par.ss.m.xu0));
cx=pipe_grad_base_s(rr,par.ss.m);
for j=1:length(rr)
    vt=zeros(size(rr)); vt(j)=k;
    cxt(j)=(pipe_obj_base_s(rr+vt,par.ss.m)-pipe_obj_base_s(rr,par.ss.m))/k;
end
norm(cx-cxt')
plot(abs([cx cxt']))
%plot(imag([cx cxt']))
%plot((cx-cxt')./(cxt+(cxt==0))')

%% test Jacobian
k=.000000001;
rr=rand(size(par.ss.m.xu0));
J1=pipe_jacobian_base(rr,par.ss.m);
J1t=sparse([]);
for j=1:length(rr)
    vt=zeros(size(rr)); vt(j)=k;
    Cip=pipe_constraints_base(rr+vt,par.ss.m); Ci=pipe_constraints_base(rr,par.ss.m);
    J1t(:,j)=sparse((Cip-Ci)/k);
    if(mod(j,100)==0) disp(num2str(j)), end
end

[norm(full(J1-J1t)) max(max(abs(J1-J1t)))]







%% Test DAE jac

[sim,part,x0]=test_jac_dae_setup(par.sim);



%% test Jacobian
k=.0000001; t=rand*part.Tgridnd(end);
rr=rand(size(x0));
drr=rand(size(x0));
[J1,J2]=pipe_net_f_imp_jac(t,rr,drr,part);
J1t=sparse([]);
for j=1:length(rr)
    vt=zeros(size(rr)); vt(j)=k;
    Cip=pipe_net_f_imp(t,rr+vt,drr,part); Ci=pipe_net_f_imp(t,rr,drr,part);
    J1t(:,j)=sparse((Cip-Ci)/k);
    if(mod(j,100)==0) disp(num2str(j)), end
end

[norm(full(J1-J1t)) max(max(abs(J1-J1t)))]

%%
%for j=1:length(xu0)
N1=par.ss.m.N1;M=par.ss.m.M;C=par.ss.m.C;FN=par.ss.m.FN;NE=par.ss.m.NE;PN=par.ss.m.PN; GN=par.ss.m.GN;
pindp=[1:FN]; pindq=[FN+1:M];
%for j=1:(N1*M)
for j=FN+1:M
%for j=N1*M+1:N1*(M+C)
%   plot([J1(j,:)-J1t(j,:)]), title(num2str(j)), pause
%   plot([J1(pindq,j)-J1t(pindq,j)]), title(num2str(j)), pause
    plot([J1(pindp,j)-J1t(pindp,j)]), title(num2str(j)), pause
%   plot((J1(pindq,j)-J1t(pindq,j))), title(num2str(j)), pause
%plot((J1(:,j)-J1t(:,j))./(J1t(:,j)+(J1t(:,j)==0))), title(num2str(j)), pause
end