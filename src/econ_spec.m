function [par]=econ_spec(par,mfolder)
%Anatoly Zlotnik, updated September 2017

mmscfd_to_kgps=par.c.mmscfd_to_kgps;
psi_to_pascal=par.c.psi_to_pascal;

if(par.intervals==0)
    %load input data on time series points
    xd=double(dlmread([mfolder '\input_ts_tpts.csv']));
    qbar=double(dlmread([mfolder '\input_ts_qbar.csv']));
    gbar=double(dlmread([mfolder '\input_ts_gbar.csv']));
    dmax_0=double(dlmread([mfolder '\input_ts_dmax.csv']));
    smax_0=double(dlmread([mfolder '\input_ts_smax.csv']));
    cd_0=double(dlmread([mfolder '\input_ts_cd.csv']));
    cs_0=double(dlmread([mfolder '\input_ts_cs.csv']));
    cslack=double(dlmread([mfolder '\input_ts_cslack.csv']));
    pslack=double(dlmread([mfolder '\input_ts_pslack.csv']));
    
    dmax_pos=find(sum(dmax_0,1)>0); %nonzero demand upper bound gnodes
    smax_pos=find(sum(smax_0,1)>0); %nonzero supply upper bound gnodes
    par.dmax_pos=dmax_pos; par.smax_pos=smax_pos;
    gs=par.n0.phys_node(smax_pos);   %gnodes with optimized supply
    gd=par.n0.phys_node(dmax_pos);   %gnodes with optimized demand 
    dmax=dmax_0(:,dmax_pos); dmin=zeros(size(dmax));
    smax=smax_0(:,smax_pos); smin=zeros(size(smax));
    cd=cd_0(:,dmax_pos); cs=cs_0(:,smax_pos);
    par.m.GN=length(gs)+length(gd); par.m.gs=gs; par.m.gd=gd;
    par.m.fn=qbar(1,:)'; par.m.pn=pslack(1,:); par.m.xd=xd'*3600;
    qbar=qbar(2:end,:); cslack=cslack(2:end,:); pslack=pslack(2:end,:);
    
    %reindexing of all gnodes to flow nodes (for integrating gbar into qbar)
    gall=par.n0.phys_node; gunique=unique(gall);
    gallind=zeros(size(gall)); guniqueind=zeros(size(gunique));
    for j=1:length(gallind)
         gallind(j)=find(par.m.fn==gall(j));
    end
    for j=1:length(guniqueind)
         guniqueind(j)=find(par.m.fn==gunique(j));
    end
    gtod=sparse(length(guniqueind),length(gall));   
    for j=1:length(guniqueind)     
        gtod(j,:)=(gallind==guniqueind(j));  
    end

    qbar(:,guniqueind)=qbar(:,guniqueind)+gbar*gtod';   %total baseline physical flow (withdrawal)
end

if(par.intervals>0)
    %load input data on time intervals
    int_qbar0=double(dlmread([mfolder '\input_int_qbar.csv']));
    int_gbar=double(dlmread([mfolder '\input_int_gbar.csv']));
    int_dmax0=double(dlmread([mfolder '\input_int_dmax.csv']));
    int_smax0=double(dlmread([mfolder '\input_int_smax.csv']));
    int_cd0=double(dlmread([mfolder '\input_int_cd.csv']));
    int_cs0=double(dlmread([mfolder '\input_int_cs.csv']));
    int_cslack=double(dlmread([mfolder '\input_int_cslack.csv']));
    int_pslack=double(dlmread([mfolder '\input_int_pslack.csv']));
    par.m.fn=int_qbar0(1,:)'; par.m.pn=int_pslack(1,:);
    
    %reindexing of all gnodes to flow nodes (for integrating gbar into qbar)
    gall=par.n0.phys_node; gunique=unique(gall);
    gallind=zeros(size(gall)); guniqueind=zeros(size(gunique));
    for j=1:length(gallind)
         gallind(j)=find(par.m.fn==gall(j));
    end
    for j=1:length(guniqueind)
         guniqueind(j)=find(par.m.fn==gunique(j));
    end
    gtod=sparse(length(guniqueind),length(gall));   
    for j=1:length(guniqueind)     
        gtod(j,:)=(gallind==guniqueind(j));  
    end
    
    %process into core input
    dmax_pos=find(sum(int_dmax0,1)>0); %nonzero demand upper bound gnodes
    smax_pos=find(sum(int_smax0,1)>0); %nonzero supply upper bound gnodes
    par.dmax_pos=dmax_pos; par.smax_pos=smax_pos;
    gs=par.n0.phys_node(smax_pos);   %gnodes with optimized supply
    gd=par.n0.phys_node(dmax_pos);   %gnodes with optimized demand
    int_qbar=int_qbar0(2:end,:); 
    int_qbar(:,guniqueind)=int_qbar(:,guniqueind)+int_gbar*gtod';   %total baseline physical flow (withdrawal)
    int_dmax=int_dmax0(:,dmax_pos);
    int_dmin=zeros(size(int_dmax));
    int_smax=int_smax0(:,smax_pos);
    int_smin=zeros(size(int_smax));
    int_cd=int_cd0(:,dmax_pos);
    int_cs=int_cs0(:,smax_pos);
    int_pslack=int_pslack(2:end,:); int_cslack=int_cslack(2:end,:); 
    par.m.gs=gs; par.m.gd=gd; par.m.GN=length(gs)+length(gd);
    xdpts_int=par.intervals;
    xd=[0:par.c.T/xdpts_int:par.c.T]'; par.m.xd=xd';
    
    %save for plotting these inputs later
    par.int_qbar=int_qbar; par.int_dmax=int_dmax; par.int_dmin=int_dmin;
    par.int_smax=int_smax; par.int_smin=int_smin; par.int_cd=int_cd;
    par.int_cs=int_cs; par.int_cslack=int_cslack; par.int_pslack=int_pslack;
    
    %transform time intervals to time point series
    In=par.intervals; Int=par.c.T/par.intervals;
    interp_tpts=reshape([[0 [1:In-1]*Int+Int/100];[[1:In-1]*Int-Int/100 par.c.T]],2*In,1);
    interp_qbar=kron(int_qbar,[1;1]);
    interp_dmax=kron(int_dmax,[1;1]);
    interp_dmin=kron(int_dmin,[1;1]);
    interp_smax=kron(int_smax,[1;1]);
    interp_smin=kron(int_smin,[1;1]);
    interp_cd=kron(int_cd,[1;1]);
    interp_cs=kron(int_cs,[1;1]);
    interp_cslack=kron(int_cslack,[1;1]);
    interp_pslack=kron(int_pslack,[1;1]);
    qbar=interp1qr(interp_tpts,interp_qbar,par.m.xd');
    dmax=interp1qr(interp_tpts,interp_dmax,par.m.xd');
    dmin=interp1qr(interp_tpts,interp_dmin,par.m.xd');
    smax=interp1qr(interp_tpts,interp_smax,par.m.xd');
    smin=interp1qr(interp_tpts,interp_smin,par.m.xd');
    cd=interp1qr(interp_tpts,interp_cd,par.m.xd');
    cs=interp1qr(interp_tpts,interp_cs,par.m.xd');
    cslack=interp1qr(interp_tpts,interp_cslack,par.m.xd');
    pslack=interp1qr(interp_tpts,interp_pslack,par.m.xd');
    
    par.int_qbar_nd=int_qbar/par.c.qsc; 
    par.int_dmax_nd=int_dmax/par.c.qsc; 
    par.int_dmin_nd=int_dmin/par.c.qsc;
    par.int_smax_nd=int_smax/par.c.qsc; 
    par.int_smin_nd=int_smin/par.c.qsc; 
    par.int_pslack_nd=int_pslack/par.c.psc;
    if(par.units==1)
        par.int_qbar_nd=par.int_qbar_nd*mmscfd_to_kgps;
        par.int_dmax_nd=par.int_dmax_nd*mmscfd_to_kgps;
        par.int_dmin_nd=par.int_dmin_nd*mmscfd_to_kgps;
        par.int_smax_nd=par.int_smax_nd*mmscfd_to_kgps;
        par.int_smin_nd=par.int_smin_nd*mmscfd_to_kgps;
        par.int_cd=int_cd/mmscfd_to_kgps;
        par.int_cs=int_cs/mmscfd_to_kgps; 
        par.int_cslack=int_cslack/mmscfd_to_kgps;
        par.int_pslack_nd=par.int_pslack_nd*psi_to_pascal;
    end
end

if(par.m.use_init_state==1)
par.m.state_init=double(dlmread([mfolder '\input_state-init.csv'])); end

%rescale and allocate input data
par.m.Yq=[sparse(qbar)'/par.c.qsc;sparse(par.n.nv-par.n0.nv,length(xd))];
par.m.Yubd=dmax'/par.c.qsc;
par.m.Ylbd=dmin'/par.c.qsc;
par.m.Yubs=smax'/par.c.qsc;
par.m.Ylbs=smin'/par.c.qsc;
par.m.Prd=cd';                       %demand price bid
par.m.Prs=cs';                       %supply price offer
par.m.Prslack=cslack';               %slack node price
par.m.Pslack=pslack'/par.c.psc;      %slack node pressure

if(par.units==1)
    par.m.Yq=par.m.Yq*mmscfd_to_kgps;
    par.m.Yubd=par.m.Yubd*mmscfd_to_kgps;
    par.m.Ylbd=par.m.Ylbd*mmscfd_to_kgps;
    par.m.Yubs=par.m.Yubs*mmscfd_to_kgps;
    par.m.Ylbs=par.m.Ylbs*mmscfd_to_kgps;
    par.m.Prd=par.m.Prd/mmscfd_to_kgps;
    par.m.Prs=par.m.Prs/mmscfd_to_kgps;
    par.m.Prslack=par.m.Prslack/mmscfd_to_kgps;
    par.m.Pslack=par.m.Pslack*psi_to_pascal;
end
    
%averages of input functions
Yfun=@(t) interp1qr(xd,par.m.Yq',t);
par.m.Yq1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Yubd',t);
par.m.Yubd1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Ylbd',t);
par.m.Ylbd1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Yubs',t);
par.m.Yubs1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Ylbs',t);
par.m.Ylbs1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Prd',t);
par.m.Prd1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Prs',t);
par.m.Prs1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Prslack',t);
par.m.Prslack1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);
Yfun=@(t) interp1qr(xd,par.m.Pslack',t);
par.m.Pslack1=quadv(@(t) Yfun(t'),0,xd(end))'/xd(end);

%reindexing of optimized gnodes to flow nodes
par.m.gall=[par.m.gd;par.m.gs];
par.m.gallsign=[ones(size(par.m.gd));-ones(size(par.m.gs))];
par.m.gunique=unique(par.m.gall);
par.m.gallind=zeros(size(par.m.gall));
par.m.guniqueind=zeros(size(par.m.gunique));
for j=1:length(par.m.gallind)
     par.m.gallind(j)=find(par.m.fn==par.m.gall(j));
end
for j=1:length(par.m.guniqueind)
     par.m.guniqueind(j)=find(par.m.fn==par.m.gunique(j));
end

%reindexing of optimized gnodes to flow nodes
par.m.gtod=sparse(length(par.m.guniqueind),par.m.GN);   
for j=1:length(par.m.guniqueind)     
    %par.m.gtod(j,:)=(par.m.gallind==par.m.guniqueind(j)); 
    par.m.gtod(j,:)=(par.m.gallind==par.m.guniqueind(j)).*par.m.gallsign;
end