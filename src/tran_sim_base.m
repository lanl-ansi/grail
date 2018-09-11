function [sim]=tran_sim_base(sim)

%If do Z
%if(sim.m.doZ==1), sim.m.Ys=p_to_rho_nd(sim.m.Ys,sim.c.b1,sim.c.b2,sim.c.psc); 
%   sim.m.Pslack1=p_to_rho_nd(sim.m.Pslack1,sim.c.b1,sim.c.b2,sim.c.psc); end

%startup simulation
t_int=[0 sim.c.T*sim.startup];
c_int=[sim.m.ccc0'; sim.m.cc0(:,1)'];
s_int=[sim.m.Pslack1'; sim.m.Ys(:,1)'];
d_int=[sim.m.Yd1';sim.m.Yd(:,1)'];
sim.m.cfun=@(t) interp1qr(t_int',c_int,t)';
sim.m.dout=@(t) interp1qr(t_int',d_int,t)';
sim.m.sfun=@(t) interp1qr(t_int',s_int,t)';
sim.m.dsfun=@(t) interp1qr(t_int',zeros(size(s_int)),t)';
sim.m.dout0=@(t) interp1qr(t_int',d_int,t)';

disp('startup simulation... ')
sim.m.Tgridnd=sim.startupgrid; sim.m.Tsc=sim.c.Tsc;
sim.m.rtol=sim.rtol0; sim.m.atol=sim.atol0;
tic,[tt0,ppnd0,qqnd0]=pipe_sim_net(sim.m); sim.init_time=toc;
disp(['finished in ' num2str(sim.init_time) ' sec.'])
sim.m.ppp0=ppnd0(end,:)'; sim.m.qqq0=qqnd0(end,:)';
cft0=sim.m.cfun(tt0)';
dft0=sim.m.dout(tt0)';

%simulation cycles
sim.m.cfun=@(t) interp1qr((-sim.m.x+1)*sim.c.T/2,sim.m.cc0',t)';
sim.m.dout=@(t) interp1qr(sim.m.xd',sim.m.Yd',t)';
sim.m.sfun=@(t) interp1qr(sim.m.xd',sim.m.Ys',t)';
xoff=mean(diff(sim.m.xd))/100;
zm=kron(ones(length(sim.m.xd)-1,1),[0;-1])*xoff;zm(end)=0;
zp=kron(ones(length(sim.m.xd)-1,1),[1;0])*xoff;zp(1)=0;
dsx=kron(sim.m.xd',[1;1]); dsx([1;length(dsx)])=[];
dsx=dsx+zm+zp;
dsy=kron(diff(sim.m.Ys')./diff(sim.m.xd'),[1;1]);
sim.m.dsfun=@(t) interp1qr(dsx,dsy,t)';

disp(['simulating first period... '])
sim.m.Tgridnd=sim.periodgrid;
sim.m.rtol=sim.rtol1; sim.m.atol=sim.atol1;
tic,[tt1,ppnd1,qqnd1]=pipe_sim_net(sim.m); sim.cycles_time=toc;
disp(['finished in ' num2str(sim.cycles_time) ' sec.'])
sim.m.ppp0=ppnd1(end,:)'; sim.m.qqq0=qqnd1(end,:)';
cft1=sim.m.cfun(mod(tt1,sim.c.T))';
dft1=sim.m.dout(mod(tt1,sim.c.T))'; 

disp(['simulating ' num2str(sim.nperiods-1) ' more periods... '])
sim.m.Tgridnd=sim.cyclesgrid;
sim.m.rtol=sim.rtol; sim.m.atol=sim.atol;
tic,[tt2,ppnd2,qqnd2]=pipe_sim_net(sim.m); sim.cycles_time=toc;
disp(['finished in ' num2str(sim.cycles_time) ' sec.'])
sim.m.ppp0=ppnd2(end,:)'; sim.m.qqq0=qqnd2(end,:)';
cft2=sim.m.cfun(mod(tt2,sim.c.T))';
dft2=sim.m.dout(mod(tt2,sim.c.T))'; 

sim.tt=[tt0(1:end-1)-sim.c.T*sim.startup;tt1(1:end-1);tt2+sim.c.T]; 
sim.ppnd=[ppnd0(1:end-1,:);ppnd1(1:end-1,:);ppnd2];
sim.qqnd=[qqnd0(1:end-1,:);qqnd1(1:end-1,:);qqnd2];
sim.cc=[cft0(1:end-1,:);cft1(1:end-1,:);cft2];
numd=sim.n0.nv-sum(sim.n0.isslack);
sim.dd=[dft0(1:end-1,1:numd);dft1(1:end-1,1:numd);dft2(:,1:numd)]*sim.c.qsc;

%process dimensional solution (metric)
sim.pp=zeros(length(sim.tt),sim.n0.nv);
%s=kron(sim.m.p_min_nd(sim.m.snodes),ones(1,length(sim.tt)));
s_nodal=sim.m.sfun(mod(sim.tt,sim.c.T));
if(sim.m.doZ==1), s_nodal=p_to_rho_nd(s_nodal,sim.c.b1,sim.c.b2,sim.c.psc); end
scf=find(ismember(sim.m.snodes,sim.m.comp_pos(sim.m.spos,1)));
s(scf,:)=sim.cc(:,sim.m.spos)'.*s_nodal(scf,:);
sim.pp(:,sim.m.snodes)=s';
sim.pp(:,sim.m.dnodes)=sim.ppnd;
sim.qq=sim.qqnd*sim.m.Xs*sim.c.qsc; 
ff=sim.qqnd*sim.c.qsc; 
for j=1:length(sim.m.dpos)
    cpj=sim.n.comp_pos(sim.m.dpos(j),1);
    sim.pp(:,cpj)=sim.pp(:,cpj).*sim.cc(:,sim.m.dpos(j));
end
%nodal pressure (before compressors)
pp_nodal=zeros(length(sim.tt),sim.n.nv);
pp_nodal(:,sim.m.snodes)=s_nodal';
pp_nodal(:,sim.m.dnodes)=sim.ppnd;

if(sim.m.doZ==1), sim.pp=rho_to_p_nd(sim.pp,sim.c.b1,sim.c.b2,sim.c.psc); 
    pp_nodal=rho_to_p_nd(pp_nodal,sim.c.b1,sim.c.b2,sim.c.psc); end


sim.pnodin=pp_nodal*sim.c.psc;   %nodal pressures (before compression)
sim.pnodout=sim.pp*sim.c.psc;     %nodal pressures (after compression)
sim.qqin=sim.qq(:,sim.n.from_flows);
sim.qqout=sim.qq(:,sim.n.to_flows);
%pipe inlet and outlet pressure (compressors only at inlets)
for j=1:sim.n0.ne
    if(sim.n0.comp_bool(j)==1)
        sim.ppin(:,j)=sim.pnodout(:,sim.n0.from_id(j));
        sim.ppout(:,j)=sim.pnodin(:,sim.n0.to_id(j));
    elseif(sim.n0.comp_bool(j)==0)
        sim.ppin(:,j)=sim.pnodin(:,sim.n0.from_id(j));
        sim.ppout(:,j)=sim.pnodin(:,sim.n0.to_id(j));
    end
end
%pipe inlet and outlet flux
sim.ffin=ff(:,sim.n.from_flows);
sim.ffout=ff(:,sim.n.to_flows);

function [p]=rho_to_p_nd(rho,b1,b2,psc)
p=(-b1+sqrt(b1^2+4*(b2*psc)*rho))/(2*b2*psc);

function [rho]=p_to_rho_nd(p,b1,b2,psc)
rho=p.*(b1+b2*psc*p);
