%Anatoly Zlotnik, December 2016

function [tr]=tran_opt_base(tr)

if(tr.Nvec(end)>1)
    if(tr.m.extension==0)   %no periodic extension over extended time horizon
        tr.dout=@(t) interp1qr(tr.m.xd',tr.m.Yq',t)';        %baseline flow function
        tr.dlbout=@(t) [interp1qr(tr.m.xd',tr.m.Ylbd',t)'; interp1qr(tr.m.xd',tr.m.Ylbs',t)'];     %lower bound on flow change
        tr.dubout=@(t) [interp1qr(tr.m.xd',tr.m.Yubd',t)'; interp1qr(tr.m.xd',tr.m.Yubs',t)'];     %upper bound on flow change demand
        tr.prdout=@(t) [interp1qr(tr.m.xd',tr.m.Prd',t)'; -interp1qr(tr.m.xd',tr.m.Prs',t)'];     %price bids function of buyers demand
        tr.prslout=@(t) interp1qr(tr.m.xd',tr.m.Prslack',t)';       %slack node price
        tr.pslout=@(t) interp1qr(tr.m.xd',tr.m.Pslack',t)';       %slack node pressure
    elseif(tr.m.extension>0) %periodic extension over extended time horizon
        tr.Nvec_old=tr.Nvec;
        tr.Nvec=[2.^([2:log2((tr.optintervals+tr.m.extension-1)/2)]) tr.optintervals+tr.m.extension-1];
        %tr.Nvec=tr.optintervals+tr.m.extension-1;   %temp check
        tr.m.Ts_old=tr.m.Ts;
        tr.m.Ts=(1+tr.m.extension/tr.optintervals)*tr.m.Ts;
        %tr.c.Tsc=(1+tr.m.extension/tr.optintervals)*tr.c.Tsc;
        ext_domain=[tr.m.xd'; (1+tr.m.extension/tr.optintervals)*tr.c.T]; 
        tr.dout=@(t) interp1qr(ext_domain,[tr.m.Yq'; tr.m.Yq(:,1)'],t)';        %baseline flow function
        tr.dlbout=@(t) [interp1qr(ext_domain,[tr.m.Ylbd'; tr.m.Ylbd(:,1)'],t)'; ...
            interp1qr(ext_domain,[tr.m.Ylbs'; tr.m.Ylbs(:,1)'],t)'];     %lower bound on flow change
        tr.dubout=@(t) [interp1qr(ext_domain,[tr.m.Yubd'; tr.m.Yubd(:,1)'],t)'; ...
            interp1qr(ext_domain,[tr.m.Yubs'; tr.m.Yubs(:,1)'],t)'];     %upper bound on flow change demand
        tr.prdout=@(t) [interp1qr(ext_domain,[tr.m.Prd'; tr.m.Prd(:,1)'],t)'; ...
            -interp1qr(ext_domain,[tr.m.Prs'; tr.m.Prs(:,1)'],t)'];     %price bids function of buyers demand
        tr.prslout=@(t) interp1qr(ext_domain,[tr.m.Prslack'; tr.m.Prslack(:,1)'],t)';       %slack node price
        tr.pslout=@(t) interp1qr(ext_domain,[tr.m.Pslack'; tr.m.Pslack(:,1)'],t)';       %slack node pressure
    end
elseif(tr.Nvec(end)==0)
     tr.Nvec_old=tr.Nvec;tr.m.Ts_old=tr.m.Ts;
end

%load sizes
FN=tr.m.FN;    % # of flow nodes
PN=tr.m.PN;    % # of pressure nodes
GN=tr.m.GN;    % # of gNodes
NE=tr.m.NE;    % # of edges
M=tr.m.M;      % # of state variables - edges and nodes
C=tr.m.C;      % # of compressors

if(tr.m.use_init_state==1)
    st_ind=0;
    pp_init=tr.m.state_init(st_ind+1:st_ind+FN); st_ind=st_ind+FN;
    qq_init=tr.m.state_init(st_ind+1:st_ind+NE); st_ind=st_ind+NE;
    cc_init=tr.m.state_init(st_ind+1:st_ind+C); st_ind=st_ind+C;
    fd_init=tr.m.state_init(st_ind+1:st_ind+GN); st_ind=st_ind+GN;
    d_init=tr.m.state_init(st_ind+1:st_ind+FN); st_ind=st_ind+FN;
    dlb_init=tr.m.state_init(st_ind+1:st_ind+GN); st_ind=st_ind+GN;
    dub_init=tr.m.state_init(st_ind+1:st_ind+GN); st_ind=st_ind+GN;
    prd_init=tr.m.state_init(st_ind+1:st_ind+GN); st_ind=st_ind+GN;
    prslack_init=tr.m.state_init(st_ind+1:st_ind+PN); st_ind=st_ind+PN;
    pslack_init=tr.m.state_init(st_ind+1:st_ind+PN); st_ind=st_ind+PN;
    lmp0_init=tr.m.state_init(st_ind+1:st_ind+FN); st_ind=st_ind+FN;
    mult0_pmax_init=tr.m.state_init(st_ind+1:st_ind+C); st_ind=st_ind+C;
    mult0_cmax_init=tr.m.state_init(st_ind+1:st_ind+C); st_ind=st_ind+C;
end

tic,
for v=1:length(tr.Nvec)        %iterate optimizations for increasing discretization



% Low-order derivative scheme
tr.m.N=tr.Nvec(v);                   % time discretization points
N=tr.m.N;       % # of time intervals
N1=tr.m.N+1;    % # of time collocation points

if(tr.m.N>1)
[tr.m.x,tr.m.D]=foDc(tr.m.N);       % compute first order collocation nodes x and differentiation matrix D
[tr.m.x,tr.m.w]=fonodes(tr.m.N);    % compute first order quadrature weights w
tr.m.t=(-tr.m.x+1)*tr.m.Ts/2;       % rescale collocation points on [-1,1] to time interval [0,T]
tr.m.N1=tr.m.N+1; tr.m.D=sparse(tr.m.D); 
tr.m.tk=tr.m.t(1:tr.m.N1);
tr.m.xk=tr.m.x(1:tr.m.N1);
elseif(tr.m.N==0)
    tr.m.x=0;tr.m.D=0;tr.m.w=1;tr.m.t=0;tr.m.N1=1;tr.m.tk=0;tr.m.xk=0;
end
% Initialize Optimization

%demands and prices on collocation nodes
if(tr.m.N==0)
    tr.m.d=tr.m.Yq1;
    tr.m.dlb=[tr.m.Ylbd1; tr.m.Ylbs1];
    tr.m.dub=[tr.m.Yubd1; tr.m.Yubs1];
    tr.m.prd=[tr.m.Prd1; -tr.m.Prs1];
    tr.m.prslack=tr.m.Prslack1;
    tr.m.pslack=tr.m.Pslack1;
    if(tr.m.doZ==1), tr.m.pslack=p_to_rho_nd(tr.m.Pslack1,tr.c.b1,tr.c.b2,tr.c.psc); end
else
tr.m.d=sparse(tr.dout(tr.m.tk*tr.c.Tsc));
tr.m.dlb=tr.dlbout(tr.m.tk*tr.c.Tsc);
tr.m.dub=tr.dubout(tr.m.tk*tr.c.Tsc);
tr.m.prd=tr.prdout(tr.m.tk*tr.c.Tsc);
tr.m.prslack=tr.prslout(tr.m.tk*tr.c.Tsc);
tr.m.pslack=tr.pslout(tr.m.tk*tr.c.Tsc);
if(tr.m.doZ==1), tr.m.pslack=p_to_rho_nd(tr.m.pslack,tr.c.b1,tr.c.b2,tr.c.psc); end
if(tr.m.use_init_state==1)
    tr.m.d(:,1)=d_init;
    tr.m.dlb(:,1)=dlb_init;
    tr.m.dub(:,1)=dub_init;
    tr.m.prd(:,1)=prd_init;
    tr.m.prslack(:,1)=prslack_init;
    tr.m.pslack(:,1)=pslack_init;
end
end

% initial guess
%variables e.g. x_{ij} have space in i (rows) and time in j (cols)
%stack columns x_i in a single vector
if(v==1)
    if(tr.m.use_init_state==0)
        pp0=kron(ones(N1,1),tr.m.p_min_nd(tr.n.nonslack_nodes));   %initialize density
        qq0=kron(zeros(N1,1),ones(NE,1));                %initialize flux
        cc0=kron(ones(N1,1),tr.n.c_max);                %initialize compressions
        dd0=zeros(N1*GN,1);                        %initialize variable offtakes
    elseif(tr.m.use_init_state==1)
        pp0=kron(ones(N,1),tr.m.p_min_nd(tr.n.nonslack_nodes));   %initialize density
        qq0=kron(zeros(N,1),ones(NE,1));                %initialize flux
        cc0=kron(ones(N,1),tr.n.c_max);                %initialize compressions
        dd0=zeros(N*GN,1);                        %initialize variable offtakes
    end
    xu0=[pp0;qq0;cc0;dd0];                               %stack it all in a vector
    xr=xu0+rand(size(xu0))/1000;                     %perturb (used to get Jacobian structure)
end

%define constraints 

%discharge pressure setpoints

%linear equality constraints
tr.m.Aeq=[];
tr.m.Beq=[];

%     tr.m.Beq=[pp_init;qq_init;cc_init;dd_init];
%     tr.m.Aeq=sparse(FN+NE+C+GN,(FN+NE+C+GN)*N1);
%     tr.m.Aeq(1:FN,1:FN)=sparse(eye(FN));
%     tr.m.Aeq(FN+1:FN+NE,N1*FN+1:N1*FN+NE)=sparse(eye(NE));
%     tr.m.Aeq(FN+NE+1:FN+NE+C,N1*(FN+NE)+1:N1*(FN+NE)+C)=sparse(eye(C));
%     tr.m.Aeq(FN+NE+C+1:FN+NE+C+GN,N1*(FN+NE+C)+1:N1*(FN+NE+C)+GN)=sparse(eye(GN));

%links from supplier slack nodes
slinks=tr.m.comp_pos(tr.m.spos,2);
clinks=tr.m.comp_pos(:,2);

%relax bounds on optimized flows at gnodes so that interval-averaged bounds dominate
bnd_rel=2;

if(tr.m.use_init_state==0)  %variables on time points 0 to N
    %lower bounds on variables
    lbp=kron(ones(N1,1),tr.m.p_min_nd(tr.n.nonslack_nodes));
    if(tr.m.doZ==1), 
    lbp=kron(ones(N1,1),p_to_rho_nd(tr.m.p_min_nd(tr.n.nonslack_nodes),tr.c.b1,tr.c.b2,tr.c.psc));      %on pressure
    end
    %lbq=kron(ones(N1,1),-10*max(abs(tr.n.q_max))*ones(NE,1)/tr.c.qsc); %on flow
    lbq_val=-10*max(abs(tr.n.q_max))/tr.c.qsc*ones(NE,1);
    lbq_val(clinks)=tr.m.flow_min_nd./tr.m.xs(clinks);
    lbq_val(slinks)=zeros(length(slinks),1);
    lbq=kron(ones(N1,1),lbq_val);                                       %on flow
    lbc=kron(ones(N1,1),tr.n.c_min);                                    %on compressor ratio
    lbd=reshape(tr.m.dlb,GN*N1,1);                                      %on flexible demands
    if(tr.intervals>1 && tr.m.N>0), lbd=(1/bnd_rel)*lbd.*(lbd>0)+bnd_rel*lbd.*(lbd<=0); end
    ipopt_options.lb = [lbp;lbq;lbc;lbd];                               %lower bound on all variables

    %upper bounds on variables
    ubp=kron(ones(N1,1),tr.m.p_max_nd(tr.n.nonslack_nodes));
    if(tr.m.doZ==1), 
    ubp=kron(ones(N1,1),p_to_rho_nd(tr.m.p_max_nd(tr.n.nonslack_nodes),tr.c.b1,tr.c.b2,tr.c.psc));        %on pressure
    end
    ubq_val=10*max(abs(tr.n.q_max))/tr.c.qsc*ones(NE,1);
    ubq_val(clinks)=tr.m.flow_max_nd./tr.m.xs(clinks);
    ubq=kron(ones(N1,1),ubq_val);                                         %on flow
    ubc=kron(ones(N1,1),tr.n.c_max);                                      %on compressor ratio
    ubd=reshape(tr.m.dub,GN*N1,1);                                        %on flexible demands
    if(tr.intervals>1 && tr.m.N>0), ubd=bnd_rel*ubd.*(ubd>0)+(1/bnd_rel)*ubd.*(ubd<0); end
    ipopt_options.ub = [ubp;ubq;ubc;ubd];                                 %upper bound on all variables
    
    %bounds on optimized flows averaged over time intervals
    if(tr.intervals>1 && v==length(tr.Nvec) && tr.m.N>0)
        tr.m.int_flow_const=1;
        Aineq_0=(1/2)*[sparse([1:GN*N1]',[1:GN*N1]',ones(GN*N1,1)) zeros(GN*N1,GN)]...
            +(1/2)*[zeros(GN*N1,GN) sparse([1:GN*N1]',[1:GN*N1]',ones(GN*N1,1))];
        Aineq_0(:,1:GN)=Aineq_0(:,1:GN)+Aineq_0(:,GN*N1+1:GN*N1+GN);Aineq_0(:,GN*N1+1:GN*N1+GN)=[];
        tr.m.Aineq=[sparse([],[],[],GN*N1,N1*(M+C)) Aineq_0];
        if(tr.m.extension==0)
            tr.m.bineqlb=reshape([tr.int_dmin_nd';-tr.int_smax_nd'],GN*N1,1);
            tr.m.binequb=reshape([tr.int_dmax_nd';-tr.int_smin_nd'],GN*N1,1);
        elseif(tr.m.extension>0)
            lbext=(tr.m.dlb(:,tr.intervals+1:N1)+[tr.m.dlb(:,tr.intervals+2:N1) tr.m.dlb(:,1)])/2; 
            ubext=(tr.m.dub(:,tr.intervals+1:N1)+[tr.m.dub(:,tr.intervals+2:N1) tr.m.dub(:,1)])/2; 
            tr.m.bineqlb=reshape([[tr.int_dmin_nd';-tr.int_smax_nd'] lbext],GN*N1,1);
            tr.m.binequb=reshape([[tr.int_dmax_nd';-tr.int_smin_nd'] ubext],GN*N1,1);
        end
    else tr.m.int_flow_const=0; end

    %lower bounds on constraints
    tr.m.hplsc=1;
    dischpliml=-kron(ones(N1,1),tr.m.p_max_nd(tr.n.comp_pos(:,1)));
    if(tr.m.doZ==1), 
    dischpliml=-kron(ones(N1,1),p_to_rho_nd(tr.m.p_max_nd(tr.n.comp_pos(:,1)),tr.c.b1,tr.c.b2,tr.c.psc));  %on discharge
    end
    hpliml=-kron(ones(N1,1),tr.m.boost_pow_max_nd)*tr.m.hplsc;       %on power
    ipopt_options.cl = [tr.m.Beq;zeros(M*N1,1);dischpliml;hpliml];   %lower bounds on the constraint functions.
    if(tr.intervals>1 && v==length(tr.Nvec) && tr.m.N>0)
        ipopt_options.cl=[ipopt_options.cl;tr.m.bineqlb]; end

    %upper bounds on constraints
    dischplimu=zeros(C*N1,1);                                        %on discharge
    hplimu=zeros(C*N1,1);                                              %on power
    ipopt_options.cu = [tr.m.Beq;zeros(M*N1,1);dischplimu;hplimu];   %upper bounds on the constraint functions.
    if(tr.intervals>1 && v==length(tr.Nvec) && tr.m.N>0)
        ipopt_options.cu=[ipopt_options.cu;tr.m.binequb]; end

    if(v==1)
    %xu0=full(min(max(2*rand(size(xu0)),ipopt_options.lb),ipopt_options.ub));
    xu0=full(min(max(xu0,ipopt_options.lb),ipopt_options.ub));
    end
elseif(tr.m.use_init_state==1) %variables on time points 1 to N
    %lower bounds on variables
    lbp=kron(ones(N,1),tr.m.p_min_nd(tr.n.nonslack_nodes));
    if(tr.m.doZ==1), 
    lbp=kron(ones(N,1),p_to_rho_nd(tr.m.p_min_nd(tr.n.nonslack_nodes),tr.c.b1,tr.c.b2,tr.c.psc));      %on pressure
    end
    %lbq=kron(ones(N1,1),-10*max(abs(tr.n.q_max))*ones(NE,1)/tr.c.qsc); %on flow
    lbq_val=-10*max(abs(tr.n.q_max))/tr.c.qsc*ones(NE,1);
    lbq_val(clinks)=tr.m.flow_min_nd./tr.m.xs(clinks);
    lbq_val(slinks)=zeros(length(slinks),1);
    lbq=kron(ones(N,1),lbq_val);                                       %on flow
    lbc=kron(ones(N,1),tr.n.c_min);                                    %on compressor ratio
    lbd=reshape(tr.m.dlb(:,2:N1),GN*N,1);                                      %on flexible demands
    if(tr.intervals>1 && tr.m.N>0), lbd=(1/bnd_rel)*lbd.*(lbd>0)+bnd_rel*lbd.*(lbd<=0); end
    ipopt_options.lb = [lbp;lbq;lbc;lbd];                               %lower bound on all variables

    %upper bounds on variables
    ubp=kron(ones(N,1),tr.m.p_max_nd(tr.n.nonslack_nodes));
    if(tr.m.doZ==1), 
    ubp=kron(ones(N,1),p_to_rho_nd(tr.m.p_max_nd(tr.n.nonslack_nodes),tr.c.b1,tr.c.b2,tr.c.psc));        %on pressure
    end
    %ubq=kron(ones(N,1),10*max(abs(tr.n.q_max))*ones(NE,1)/tr.c.qsc);     %on flow
    ubq_val=10*max(abs(tr.n.q_max))/tr.c.qsc*ones(NE,1);
    ubq_val(clinks)=tr.m.flow_max_nd./tr.m.xs(clinks);
    ubq=kron(ones(N,1),ubq_val);
    ubc=kron(ones(N,1),tr.n.c_max);                                      %on compressor ratio
    ubd=reshape(tr.m.dub(:,2:N1),GN*N,1);                                        %on flexible demands
    if(tr.intervals>1 && tr.m.N>0), ubd=bnd_rel*ubd.*(ubd>0)+(1/bnd_rel)*ubd.*(ubd<0); end
    ipopt_options.ub = [ubp;ubq;ubc;ubd];                                 %upper bound on all variables

    %bounds on optimized flows averaged over time intervals
    if(tr.intervals>1 && v==length(tr.Nvec) && tr.m.N>0)
        tr.m.int_flow_const=1;
        Aineq_0=(1/2)*[sparse([1:GN*N1]',[1:GN*N1]',ones(GN*N1,1)) zeros(GN*N1,GN)]...
            +(1/2)*[zeros(GN*N1,GN) sparse([1:GN*N1]',[1:GN*N1]',ones(GN*N1,1))];
        Aineq_0(:,1:GN)=Aineq_0(:,1:GN)+Aineq_0(:,GN*N1+1:GN*N1+GN);Aineq_0(:,GN*N1+1:GN*N1+GN)=[];
        tr.m.Aineq=[sparse([],[],[],GN*N1,N1*(M+C)) Aineq_0];
        if(tr.m.extension==0)
            tr.m.bineqlb=reshape([tr.int_dmin_nd';-tr.int_smax_nd'],GN*N1,1);
            tr.m.binequb=reshape([tr.int_dmax_nd';-tr.int_smin_nd'],GN*N1,1);
        elseif(tr.m.extension>0)
            lbext=(tr.m.dlb(:,tr.intervals+1:N1)+[tr.m.dlb(:,tr.intervals+2:N1) tr.m.dlb(:,1)])/2; 
            ubext=(tr.m.dub(:,tr.intervals+1:N1)+[tr.m.dub(:,tr.intervals+2:N1) tr.m.dub(:,1)])/2; 
            tr.m.bineqlb=reshape([[tr.int_dmin_nd';-tr.int_smax_nd'] lbext],GN*N1,1);
            tr.m.binequb=reshape([[tr.int_dmax_nd';-tr.int_smin_nd'] ubext],GN*N1,1);
        end
    else tr.m.int_flow_const=0; end
    
    %lower bounds on constraints
    tr.m.hplsc=1;
    dischpliml=-kron(ones(N,1),tr.m.p_max_nd(tr.n.comp_pos(:,1)));
    if(tr.m.doZ==1), 
    dischpliml=-kron(ones(N,1),p_to_rho_nd(tr.m.p_max_nd(tr.n.comp_pos(:,1)),tr.c.b1,tr.c.b2,tr.c.psc));  %on discharge
    end
    hpliml=-kron(ones(N,1),tr.m.boost_pow_max_nd)*tr.m.hplsc;       %on power
    ipopt_options.cl = [tr.m.Beq;zeros(M*N,1);dischpliml;hpliml];   %lower bounds on the constraint functions.

    %upper bounds on constraints
    dischplimu=zeros(C*N,1);                                        %on discharge
    hplimu=zeros(C*N,1);                                              %on power
    ipopt_options.cu = [tr.m.Beq;zeros(M*N,1);dischplimu;hplimu];   %upper bounds on the constraint functions.

    if(v==1)
    %xu0=full(min(max(2*rand(size(xu0)),ipopt_options.lb),ipopt_options.ub));
    xu0=full(min(max(xu0,ipopt_options.lb),ipopt_options.ub));
    end
end


% resample xf

%this is used to resample the solution when discretization is increased
if(v>1)
    N1s=tr.m.N1s;
    Ns=tr.m.Ns;
    if(tr.m.use_init_state==0)
        p=reshape(xf(1:N1s*FN),FN,N1s); 
        q=reshape(xf(N1s*FN+1:N1s*M),NE,N1s);
        comps=reshape(xf(N1s*M+1:N1s*M+N1s*C),C,N1s);
        fd=reshape(xf(N1s*(M+C)+1:N1s*(M+C+GN)),GN,N1s);
        pn=interp1(tr.m.xks,p',tr.m.xk,'linear','extrap'); 
        qn=interp1(tr.m.xks,q',tr.m.xk,'linear','extrap'); 
        cn=interp1(tr.m.xks,comps',tr.m.xk,'linear','extrap');
        fdn=interp1(tr.m.xks,fd',tr.m.xk,'linear','extrap');
        xu0=[reshape(pn',FN*N1,1);reshape(qn',NE*N1,1);reshape(cn',C*N1,1);reshape(fdn',GN*N1,1)];
    elseif(tr.m.use_init_state==1)
        p=[pp_init reshape(xf(1:Ns*FN),FN,Ns)]; 
        q=[qq_init reshape(xf(Ns*FN+1:Ns*M),NE,Ns)];
        comps=[cc_init reshape(xf(Ns*M+1:Ns*M+Ns*C),C,Ns)];
        fd=[fd_init reshape(xf(Ns*(M+C)+1:Ns*(M+C+GN)),GN,Ns)];
        pn=interp1(tr.m.xks,p',tr.m.xk,'linear','extrap'); 
        qn=interp1(tr.m.xks,q',tr.m.xk,'linear','extrap'); 
        cn=interp1(tr.m.xks,comps',tr.m.xk,'linear','extrap');
        fdn=interp1(tr.m.xks,fd',tr.m.xk,'linear','extrap');
        xu0=[reshape(pn(2:N1,:)',FN*N,1);reshape(qn(2:N1,:)',NE*N,1);...
        reshape(cn(2:N1,:)',C*N,1);reshape(fdn(2:N1,:)',GN*N,1)];
    end
    xu0=min(max(xu0,ipopt_options.lb),ipopt_options.ub);    %make sure variables are in feasible set
    xr=xu0+rand(size(xu0))/1000;   
end
 
% IPOPT parameters

%scaling parameters
tr.m.smsc=tr.m.C*tr.m.Ts/2*max(tr.Nvec)/tr.m.cdw;  %weight cost on derivative of compressions in final step
tr.m.smsd=tr.m.GN*tr.m.Ts/2*max(tr.Nvec)/tr.m.ddw;  %weight cost on derivative of flexible demands in final step
tr.m.objsc=tr.m.Ts/2/tr.m.odw;                 %this scales the objective for ipopt

%Jacobian sparsity pattern
JS=sparse(abs(sign(pipe_jacobian_base(xr,tr.m))));    
sum(sum(abs(JS)>0))/prod(size(JS))

%Set the IPOPT options.
ipopt_options.ipopt.mu_strategy = 'adaptive';
ipopt_options.ipopt.hessian_approximation = 'limited-memory';
%options.ipopt.limited_memory_update_type = 'BFGS';
ipopt_options.ipopt.limited_memory_update_type = 'sr1';
ipopt_options.ipopt.max_iter = tr.m.maxiter;
ipopt_options.ipopt.print_level=5;
ipopt_tol=tr.m.opt_tol; %if(v==length(Nvec)) ipopt_tol=1e-4; end 
ipopt_options.ipopt.tol = ipopt_tol;
ipopt_options.ipopt.constr_viol_tol = ipopt_tol/10;
ipopt_options.ipopt.acceptable_tol= ipopt_tol;
ipopt_options.ipopt.output_file=tr.output_file;
  
% The callback functions.
if(v<length(tr.Nvec))
    ipopt_funcs.objective         = @(xu) pipe_obj_base(xu,tr.m);
    ipopt_funcs.gradient          = @(xu) pipe_grad_base(xu,tr.m);
else
   ipopt_funcs.objective         = @(xu) pipe_obj_base_s(xu,tr.m);
   ipopt_funcs.gradient          = @(xu) pipe_grad_base_s(xu,tr.m);
end
ipopt_funcs.constraints       = @(xu) pipe_constraints_base(xu,tr.m);
ipopt_funcs.jacobian          = @(xu) pipe_jacobian_base(xu,tr.m);
ipopt_funcs.jacobianstructure = @() JS;
ipopt_funcs.hessian           = @() 1;
ipopt_funcs.hessianstructure  = @() 1;

 
% run optimization

disp(['Solving with ' num2str(tr.Nvec(v)) ' points...'])
  % Run IPOPT several times (can't exit from mex file, so limit iterations)
  [xf,ip_info] = ipopt(xu0,ipopt_funcs,ipopt_options);
 if(ip_info.iter>ipopt_options.ipopt.max_iter-2)
     disp(['Solving again with ' num2str(tr.Nvec(v)) ' points...'])
   [xf,ip_info] = ipopt(xf,ipopt_funcs,ipopt_options);
end
 
%% Save old stuff
tr.m.xks=tr.m.xk;       %old LGL time grid for resampling
tr.m.N1s=tr.m.N1;    	%old time points
tr.m.Ns=tr.m.N;         %old time points
tr.xf=xf;

end
ipopt_net_time=toc, tr.m.ipopt_net_time=ipopt_net_time;

% save output

xf=tr.xf;
tr.objval=pipe_obj_base(xf,tr.m);
tr.resid=pipe_constraints_base(xf,tr.m);
tr.ip_info=ip_info;

if(tr.m.use_init_state==0)
    %restriction to original horizon if extension used
    tr.pp0=[reshape(xf(1:N1*FN),FN,N1) xf(1:FN)]; 
    tr.qq0=[reshape(xf(N1*FN+1:N1*M),NE,N1) xf(N1*FN+1:N1*FN+NE)];
    tr.cc0=[reshape(xf(N1*M+1:N1*M+N1*C),C,N1) xf(N1*M+1:N1*M+C)];
    tr.fd0=[reshape(xf(N1*(M+C)+1:N1*(M+C+GN)),GN,N1) xf(N1*(M+C)+1:N1*(M+C)+GN)];
    tr.lmp0=-[reshape(ip_info.lambda(1:FN*N1),FN,N1) ip_info.lambda(1:FN)]/tr.m.odw;
    tr.mult0_pmax=[reshape(ip_info.lambda(N1*M+1:N1*(M+C)),C,N1) ip_info.lambda(N1*M+1:N1*M+C)]/tr.m.odw;
    tr.mult0_cmax=[reshape(ip_info.lambda(N1*(M+C)+1:N1*(M+2*C)),C,N1) ip_info.lambda(N1*(M+C)+1:N1*(M+C)+C)]/tr.m.odw;
elseif(tr.m.use_init_state==1)
    %restriction to original horizon if extension used
    tr.pp0=[pp_init reshape(xf(1:N*FN),FN,N) xf(1:FN) pp_init]; 
    tr.qq0=[qq_init reshape(xf(N*FN+1:N*M),NE,N) xf(N*FN+1:N*FN+NE) qq_init];
    tr.cc0=[cc_init reshape(xf(N*M+1:N*M+N*C),C,N) xf(N*M+1:N*M+C) cc_init];
    tr.fd0=[fd_init reshape(xf(N*(M+C)+1:N*(M+C+GN)),GN,N) xf(N*(M+C)+1:N*(M+C)+GN) fd_init];
    tr.lmp0=[lmp0_init -reshape(ip_info.lambda(1:FN*N),FN,N)/tr.m.odw lmp0_init];
    tr.mult0_pmax=[mult0_pmax_init reshape(ip_info.lambda(N*M+1:N*(M+C)),C,N)/tr.m.odw mult0_pmax_init];
    tr.mult0_cmax=[mult0_cmax_init reshape(ip_info.lambda(N*(M+C)+1:N*(M+2*C)),C,N)/tr.m.odw mult0_cmax_init];
end
if(tr.m.extension>0)
    %re-scaling for LMPs to eliminate extension effect (?)
    %tr.m1=tr.m;tr.m1.w(tr.Nvec_old(end)+1:end)=0;
    %tr.lmpsc=pipe_obj_base(xf,tr.m)/pipe_obj_base(xf,tr.m1);
    tr.lmpsc=tr.Nvec(end)/tr.Nvec_old(end);
    
    %restriction onto interval of interest
    tr.m.Ts=tr.m.Ts_old; tr.Nvec=tr.Nvec_old;
    tr.m.N1=tr.Nvec_old(end)+1; N1=tr.m.N1;
    %restrict final inputs
    tr.m.d=tr.m.d(:,1:N1+1);
    tr.m.dlb=tr.m.dlb(:,1:N1+1);
    tr.m.dub=tr.m.dub(:,1:N1+1);
    tr.m.prd=tr.m.prd(:,1:N1+1);
    tr.m.prslack=tr.m.prslack(:,1:N1+1);
    tr.m.pslack=tr.m.pslack(:,1:N1+1);
    %restrict final solutions
    tr.pp0=tr.pp0(:,1:N1+1);
    tr.qq0=tr.qq0(:,1:N1+1);
    tr.cc0=tr.cc0(:,1:N1+1);
    tr.fd0=tr.fd0(:,1:N1+1);
    tr.lmp0=tr.lmp0(:,1:N1+1);
    tr.mult0_pmax=tr.mult0_pmax(:,1:N1+1);
    tr.mult0_cmax=tr.mult0_cmax(:,1:N1+1);
    tr.ip_info=ip_info;
    
    %apply re-scaling
    %     tr.lmp0=tr.lmp0*tr.lmpsc;
    %     tr.mult0_pmax=tr.mult0_pmax*tr.lmpsc;
    %     tr.mult0_cmax=tr.mult0_cmax*tr.lmpsc;
elseif(tr.m.extension==0)
    tr.m.d=[tr.m.d tr.m.d(:,1)];
    tr.m.dlb=[tr.m.dlb tr.m.dlb(:,1)];
    tr.m.dub=[tr.m.dub tr.m.dub(:,1)];
    tr.m.prd=[tr.m.prd tr.m.prd(:,1)];
    tr.m.prslack=[tr.m.prslack tr.m.prslack(:,1)];
    tr.m.pslack=[tr.m.pslack tr.m.pslack(:,1)];
end
tr.tt0=[0:tr.m.N1]'*tr.c.T/tr.m.N1;
    
if(tr.m.save_state==1)
    tr.state_save=[...
        tr.pp0(:,tr.m.state_save_pts);...
        tr.qq0(:,tr.m.state_save_pts);...
        tr.cc0(:,tr.m.state_save_pts);...
        tr.fd0(:,tr.m.state_save_pts);...
        tr.m.d(:,tr.m.state_save_pts);...
        tr.m.dlb(:,tr.m.state_save_pts);...
        tr.m.dub(:,tr.m.state_save_pts);...
        tr.m.prd(:,tr.m.state_save_pts);...
        tr.m.prslack(:,tr.m.state_save_pts);...
        tr.m.pslack(:,tr.m.state_save_pts);...
        tr.lmp0(:,tr.m.state_save_pts);... 
        tr.mult0_pmax(:,tr.m.state_save_pts);...
        tr.mult0_cmax(:,tr.m.state_save_pts)];
end


if(tr.Nvec(end)>1),
tr=rmfield(tr,{'dout','dlbout','dubout','prdout','prslout'}); end



%saved state uses density
ppp0=tr.pp0(:,1); qqq0=tr.qq0(:,1); ccc0=tr.cc0(:,1);
tr.m.ppp0=ppp0; tr.m.qqq0=qqq0; tr.m.ccc0=ccc0;




function [x,D]=foDc(N)

%x=-([0:N]'-N/2)/(N/2);
x=-([0:N+1]'-(N+1)/2)/((N+1)/2);
D=eye(N+1)-diag(ones(N,1),-1); 
D(1,N+1)=-1;
D=-D*(N/2);

function [x,w]=fonodes(N)

%x=-([0:N]'-N/2)/(N/2);
x=-([0:N+1]'-(N+1)/2)/((N+1)/2);
w=ones(N+1,1)/(N/2);

function [p]=rho_to_p_nd(rho,b1,b2,psc)
p=(-b1+sqrt(b1^2+4*b2*psc*rho))/(2*b2*psc);

function [rho]=p_to_rho_nd(p,b1,b2,psc)
rho=p.*(b1+b2*psc*p);
