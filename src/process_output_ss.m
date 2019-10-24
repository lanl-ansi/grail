function [par]=process_output_ss(par)
% Anatoly Zlotnik, February 2019

% mfolder=par.mfolder;
out=par.out;
% 
ss=par.ss; 
psi_to_pascal=ss.c.psi_to_pascal;
mpa_to_psi=1000000/psi_to_pascal;
ss.c.mpa_to_psi=mpa_to_psi;
mmscfd_to_kgps=ss.c.mmscfd_to_kgps;
hp_to_watt=745.7;
if(par.out.doZ==1), b1=ss.c.b1; b2=ss.c.b2; end

%process optimization output
%process dimensional solution (metric)
pp=zeros(length(ss.tt0),ss.n0.nv);
%s_nodal=[ss.m.pslack ss.m.pslack(:,1)];
s_nodal=ss.m.pslack;
if(par.ss.m.extension>0), s_nodal=ss.m.pslack; end
scf=find(ismember(ss.m.snodes,ss.m.comp_pos(ss.m.spos,1)));
s(scf,:)=ss.cc0(ss.m.spos,:).*s_nodal(scf,:);
pp(:,ss.m.snodes)=s';
pp(:,ss.m.dnodes)=ss.pp0';
qq=ss.qq0'*ss.m.Xs*ss.c.qsc; 
ff=ss.qq0'*ss.c.qsc; 
for j=1:length(ss.m.dpos)
    cpj=ss.n.comp_pos(ss.m.dpos(j),1);
    pp(:,cpj)=pp(:,cpj).*ss.cc0(ss.m.dpos(j),:)';
end
%nodal pressure (before compressors)
pp_nodal=zeros(length(ss.tt0),ss.n.nv);
pp_nodal(:,ss.m.snodes)=s_nodal';
pp_nodal(:,ss.m.dnodes)=ss.pp0';

if(ss.m.doZ==1), pp=rho_to_p_nd(pp,ss.c.b1,ss.c.b2,ss.c.psc);
    pp_nodal=rho_to_p_nd(pp_nodal,ss.c.b1,ss.c.b2,ss.c.psc); end

%compute mass in pipes
if(ss.m.doZ==1), p_density=p_to_rho(pp',ss.c.b1,ss.c.b2,par.ss);
else p_density=ss.c.psc*pp'/(ss.c.gasR*ss.c.gasT); end
p_comp=par.ss.cc0;
out.ss.p_mass=pipe_mass(p_density,p_comp,par.ss.m);    %all pipes
out.ss.pipe_mass_0=(par.ss.n.disc_to_edge*out.ss.p_mass)';       %original pipes

ss.pnodin=pp_nodal*ss.c.psc;   %nodal pressures (before compression)
ss.pnodout=pp*ss.c.psc;     %nodal pressures (after compression)
ss.qqin=qq(:,ss.n.from_flows);
ss.qqout=qq(:,ss.n.to_flows);
%pipe inlet and outlet pressure (compressors only at inlets)
for j=1:ss.n0.ne
    if(ss.n0.comp_bool(j)==1)
        %ss.ppin(:,j)=ss.pnodout(:,ss.n0.from_id(j));
        ss.ppin(:,j)=ss.pnodin(:,ss.n0.from_id(j));
        ss.ppout(:,j)=ss.pnodin(:,ss.n0.to_id(j));
    elseif(ss.n0.comp_bool(j)==0)
        ss.ppin(:,j)=ss.pnodin(:,ss.n0.from_id(j));
        ss.ppout(:,j)=ss.pnodin(:,ss.n0.to_id(j));
    end
end
%pipe inlet and outlet flux
ss.ffin=ff(:,ss.n.from_flows);
ss.ffout=ff(:,ss.n.to_flows);

out.ss.tt0=ss.tt0/3600;         %plotting time in hours
out.ss.qqinopt=ss.qqin;                  %flow boundary in
out.ss.qqoutopt=ss.qqout;                %flow boundary out
%if(par.out.plotnodal==1)
    out.ss.ppoptnodal=ss.pnodin(:,1:ss.n0.nv);
    out.ss.ppopt=ss.pnodout(:,1:ss.n0.nv);  %pressure (nodal)
    out.ss.qqopt=[out.ss.qqinopt out.ss.qqoutopt];  %all boundary flows
%else
    %out.ppoptall=ss.pnodout(:,1:ss.n.nv);   %pressure (all)
    %out.qqoptall=ss.qq0;                      %flows (all)
%end
    out.ss.ppinopt=ss.ppin; out.ss.ppoutopt=ss.ppout;
if(par.out.units==1), out.ss.ppopt=out.ss.ppopt/psi_to_pascal; out.ss.ppoptnodal=out.ss.ppoptnodal/psi_to_pascal; 
out.ss.ppinopt=ss.ppin/psi_to_pascal; out.ss.ppoutopt=ss.ppout/psi_to_pascal; 
out.ss.qqopt=out.ss.qqopt/mmscfd_to_kgps; out.ss.qqinopt=out.ss.qqinopt/mmscfd_to_kgps;  
out.ss.qqoutopt=out.ss.qqoutopt/mmscfd_to_kgps; out.ss.pipe_mass_0=out.ss.pipe_mass_0/mmscfd_to_kgps/86400; end

%market flow solution
ss.m.Yd=[ss.m.Yq1(1:ss.m.FN) ss.m.Yq1(1:ss.m.FN)];
ss.m.Yd(ss.m.guniqueind,:)=ss.m.Yd(ss.m.guniqueind,:)+ss.m.gtod*ss.fd0;
%ss.m.Yd=interp1qr(ss.m.xd',ss.m.Yq(1:ss.m.FN,:)',ss.tt0)';
%ss.m.Yd(ss.m.guniqueind,:)=ss.m.Yd(ss.m.guniqueind,:)+ss.m.gtod*ss.fd0;
ss.m.Ygd=ss.fd0(1:length(ss.m.gd),:);
ss.m.Ygs=-ss.fd0(length(ss.m.gd)+1:length(ss.m.gall),:);

%compressor discharge pressures
%if(par.out.dosim==1), out.csetsim=psim1(:,ss.m.comp_pos(:,1)); end
out.ss.csetopt=out.ss.ppopt(:,ss.m.comp_pos(:,1));

%process parameters (compression ratios and demands)
out.ss.cc=ss.cc0'; out.ss.td=ss.m.xd/3600; 
out.ss.dbase=ss.m.Yq(1:length(ss.m.fn),:)'*ss.c.qsc;     %base flow "q"
out.ss.gsub=ss.m.Yubs'*ss.c.qsc;   %upper bounds on sales
out.ss.gslb=ss.m.Ylbs'*ss.c.qsc;   %lower bounds on sales
out.ss.gdub=ss.m.Yubd'*ss.c.qsc;   %upper bounds on buys
out.ss.gdlb=ss.m.Ylbd'*ss.c.qsc;   %lower bounds on buys
out.ss.gdsol=ss.m.Ygd'*ss.c.qsc;   %gnode buyer solutions
out.ss.gssol=-ss.m.Ygs'*ss.c.qsc;   %gnode seller solutions
%gnode buyer and seller solutions for all original gnodes 
GN0=length(ss.n0.phys_node);    %number of original gnodes
out.ss.gdsol_all=zeros(2,GN0); out.ss.gdsol_all(:,ss.dmax_pos)=out.ss.gdsol;
out.ss.gssol_all=zeros(2,GN0); out.ss.gssol_all(:,ss.smax_pos)=out.ss.gssol;
out.ss.dgflows=full(ss.m.Yd(ss.m.guniqueind,:))'*ss.c.qsc; %flow at nodes with gnodes
slinks=ss.m.comp_pos(ss.m.spos,2); 
out.ss.supp_flow=qq(:,slinks);    %supply flow 
out.ss.nonslack_flow=full(ss.m.Yd(1:length(ss.n0.nonslack_nodes),:))*ss.c.qsc;
out.ss.flows_all=zeros(ss.n0.nv,par.ss.m.N1+1);
out.ss.flows_all(ss.n0.slack_nodes,:)=-out.ss.supp_flow;
out.ss.flows_all(ss.n0.nonslack_nodes,:)=out.ss.nonslack_flow;
out.ss.dgflows_all=out.ss.flows_all(ss.n0.phys_node,:); %flow at all original gnodes
% out.ss.supp_flow_sim=sim.qq(:,slinks);    %supply flow 
% out.ss.flows_all_sim=zeros(sim.n0.nv,length(sim.tt));
% out.ss.flows_all_sim(sim.n0.slack_nodes,:)=-out.ss.supp_flow_sim;
% out.ss.flows_all_sim(sim.n0.nonslack_nodes,:)=full(sim.m.Yd(1:length(sim.n0.nonslack_nodes),:))*sim.c.qsc;
if(par.out.units==1), out.ss.dbase=out.ss.dbase/mmscfd_to_kgps; out.ss.gsub=out.ss.gsub/mmscfd_to_kgps; 
    out.ss.gslb=out.ss.gslb/mmscfd_to_kgps; out.ss.gdub=out.ss.gdub/mmscfd_to_kgps; out.ss.gdlb=out.ss.gdlb/mmscfd_to_kgps;
    out.ss.gdsol=out.ss.gdsol/mmscfd_to_kgps; out.ss.gssol=out.ss.gssol/mmscfd_to_kgps;
    out.ss.gdsol_all=out.ss.gdsol_all/mmscfd_to_kgps; out.ss.gssol_all=out.ss.gssol_all/mmscfd_to_kgps;
    out.ss.dgflows=out.ss.dgflows/mmscfd_to_kgps; out.ss.dgflows_all=out.ss.dgflows_all/mmscfd_to_kgps; 
    out.ss.dgflows_all=out.ss.dgflows_all/mmscfd_to_kgps;  out.ss.supp_flow=out.ss.supp_flow/mmscfd_to_kgps; 
    out.ss.flows_all=out.ss.flows_all/mmscfd_to_kgps;  
    %out.ss.supp_flow_sim/mmscfd_to_kgps; %out.ss.flows_all_sim/mmscfd_to_kgps;  
end

%check nodal flow balance
out.ss.flowbal=ss.n0.Amp*out.ss.qqoutopt'+ss.n0.Amm*out.ss.qqinopt'-out.ss.flows_all;
out.ss.flowbalrel=3*out.ss.flowbal./(abs(ss.n0.Amp*out.ss.qqoutopt')+abs(ss.n0.Amm*out.ss.qqinopt')+abs(out.ss.flows_all));
out.ss.flowbalrel(mean(out.ss.flowbal')./mean(out.ss.flowbalrel')<ss.m.opt_tol,:)=0; out.ss.flowbalrel=out.ss.flowbalrel';

%out.ss.flowbals=ss.n0.Amp*out.ss.qqoutsim'+ss.n0.Amm*out.ss.qqinsim'-out.ss.flows_all;
%out.ss.flowbalsrel=3*out.ss.flowbal./(abs(ss.n0.Amp*out.ss.qqoutsim')+abs(ss.n0.Amm*out.ss.qqinsim')+abs(out.ss.flows_all));
%out.ss.flowbalsrel(mean(out.ss.flowbal')./mean(out.ss.flowbalrel')<ss.m.opt_tol,:)=0;


%compressor power
% if(par.out.dosim==1)
%     cpossim=par.ss.m.comp_pos; m=ss.m.mpow;
%     out.ss.cccom=interp1qr(out.ss.tt0,ss.cc0',out.ss.ttcom);
%     qcompsim=interp1qr(out.ss.tt,sim.qq(:,cpossim(:,2)),ttcomsim); 
%     cpow_nd=(abs(qcompsim)).*((out.ss.cccom).^(2*m)-1);
%     out.ss.cpowsim=cpow_nd.*kron(ss.m.eff',ones(size(cpow_nd,1),1))*ss.c.mmscfd_to_hp/mmscfd_to_kgps;
% end
out.ss.ccopt=ss.cc0';
cposopt=par.ss.m.comp_pos; m=ss.m.mpow;
qcompopt=qq(:,cposopt(:,2)); %cpow_nd=(abs(qcompopt)).*((ss.cc0').^(m)-1);
%out.ss.cpowopt=cpow_nd.*kron(ss.m.eff',ones(size(cpow_nd,1),1))*ss.c.mmscfd_to_hp/mmscfd_to_kgps;
out.ss.cpowopt=(abs(qcompopt)).*((ss.cc0').^(m)-1)*ss.m.Wc;   %comp power in Watts

%process locational marginal price
% out.ss.lmptr=par.ss.lmp0(par.ss.m.flexnodes,:)'/2*par.ss.m.N/par.ss.m.odw*par.ss.c.Tsc*par.ss.c.Tsc/2;
% lmpss=par.ss.lmp0(par.ss.m.flexnodes,:)'/par.ss.m.odw*par.ss.c.Tsc*par.ss.c.Tsc/2;
% out.ss.lmptr=par.ss.lmp0'/2*par.ss.m.N/par.ss.m.odw*par.ss.c.Tsc*par.ss.c.Tsc/2;
% lmpss=par.ss.lmp0'/ss.m.odw*ss.c.Tsc*ss.c.Tsc/2;
if(ss.m.N==0), trmN=2; end
if(ss.m.N>1), trmN=ss.m.N; end
out.ss.trlmp=ss.lmp0'/2*trmN;    %all lmps
out.ss.trlmpnodal=out.ss.trlmp(:,1:length(ss.m.fn));
out.ss.trlmpnodal_all=zeros(par.ss.m.N1+1,ss.n0.nv);
out.ss.trlmpnodal_all(:,ss.n0.slack_nodes)=-ss.m.prslack';
out.ss.trlmpnodal_all(:,ss.n0.nonslack_nodes)=out.ss.trlmpnodal;
out.ss.gnodelmp=out.ss.trlmpnodal_all(:,ss.n0.phys_node);
%out.ss.gnodelmp=ss.lmp0(ss.n0.phys_node,:)'/2*trmN;     %lmps at all gnodes
out.ss.dglmp=ss.lmp0(ss.m.guniqueind,:)'/2*trmN;     %lmps at market nodes
out.ss.gdlmp=ss.lmp0(ss.m.gallind(1:length(ss.m.gd)),:)'/2*trmN;     %lmps at demand gnodes
out.ss.gslmp=ss.lmp0(ss.m.gallind(length(ss.m.gd):length(ss.m.gall)),:)'/2*trmN;     %lmps at supply gnodes
% if(par.out.ss.dosim==1), out.ss.lmpss=par.ss.lmp0'; end
out.ss.Prd=ss.m.Prd; out.ss.Prs=ss.m.Prs; 
out.ss.Prslack=interp1qr(ss.m.xd',ss.m.Prslack',out.ss.tt0);  %bid and offer prices
out.ss.mult0_pmax=ss.mult0_pmax'/2*trmN*ss.c.psi_to_pascal/ss.c.psc*3600;    %output pressure marginal prices ($/Psi/hr)
out.ss.mult0_cmax=ss.mult0_cmax'/2*trmN*3.6/0.75; %compression marginal prices ($/hp)
if(par.out.units==1), out.ss.trlmp=out.ss.trlmp*mmscfd_to_kgps; 
     out.ss.trlmpnodal=out.ss.trlmpnodal*mmscfd_to_kgps; 
     out.ss.dglmp=out.ss.dglmp*mmscfd_to_kgps; out.ss.gnodelmp=out.ss.gnodelmp*mmscfd_to_kgps;
     out.ss.gdlmp=out.ss.gdlmp*mmscfd_to_kgps; out.ss.gslmp=out.ss.gslmp*mmscfd_to_kgps; 
     out.ss.Prd=out.ss.Prd*mmscfd_to_kgps; out.ss.Prs=out.ss.Prs*mmscfd_to_kgps;
     out.ss.Prslack=out.ss.Prslack*mmscfd_to_kgps;
     out.ss.cpowopt=out.ss.cpowopt/hp_to_watt;
     %if(par.out.ss.dosim==1), out.ss.lmpss=out.ss.lmpss*mmscfd_to_kgps; end
end
out.ss.lmpin=zeros(size(out.ss.ppinopt)); out.ss.lmpout=zeros(size(out.ss.ppoutopt));
for j=1:ss.n0.ne
    if(ismember(ss.n0.from_id(j),ss.m.pn))
        out.ss.lmpin(:,j)=out.ss.Prslack(:,find(ss.m.pn==ss.n0.from_id(j))); else
        out.ss.lmpin(:,j)=out.ss.trlmpnodal(:,find(ss.m.fn==ss.n0.from_id(j))); end
    if(ismember(ss.n0.to_id(j),ss.m.pn))
        out.ss.lmpout(:,j)=out.ss.Prslack(:,find(ss.m.pn==ss.n0.to_id(j))); else
        out.ss.lmpout(:,j)=out.ss.trlmpnodal(:,find(ss.m.fn==ss.n0.to_id(j))); end
end
%cmap=colormap;
%set(0,'DefaultAxesColorOrder',cmap(floor(rand(length(ss.m.gs),1)*64)+1,:))

out.ss.guniqueind=ss.m.guniqueind; out.ss.gunique=ss.m.gunique; out.ss.fn=ss.m.fn; out.ss.pn=ss.m.pn;
out.ss.n0=ss.n0; out.ss.n=ss.n; out.ss.gd=ss.m.gd; out.ss.gs=ss.m.gs; out.ss.FN=ss.m.FN; out.ss.PN=ss.m.PN;
out.ss.cn=ss.m.C; 
out.ss.mfolder=par.mfolder;

if(par.out.savecsvoutput==1)
        mfolder=par.mfolder;
%     if(par.out.intervals_out==0)
        pipe_cols=[1:out.ss.n0.ne-out.ss.n0.nc]; comp_cols=[out.ss.n0.ne-out.ss.n0.nc+1:out.ss.n0.ne];
        %dlmwrite([mfolder '\output_ss_tpts.csv'],double(out.ss.tt0),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_pipe-pressure-in.csv'],double([pipe_cols;out.ss.ppinopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_pipe-pressure-out.csv'],double([pipe_cols;out.ss.ppoutopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-pressure-in.csv'],double([1:out.ss.n0.nc;out.ss.ppinopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-pressure-out.csv'],double([1:out.ss.n0.nc;out.ss.ppoutopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_pipe-flow-in.csv'],double([pipe_cols;out.ss.qqinopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_pipe-flow-out.csv'],double([pipe_cols;out.ss.qqoutopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-flow-in.csv'],double([1:out.ss.n0.nc;out.ss.qqinopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-flow-out.csv'],double([1:out.ss.n0.nc;out.ss.qqoutopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_nodal-pressure.csv'],double([[1:out.ss.n0.nv];out.ss.ppoptnodal(1,:)]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ss_gnode-physical-withdrawals.csv'],double([out.ss.gunique';out.ss.dgflows]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_nonslack-flows.csv'],double([out.ss.fn';out.ss.flows_all(ss.n0.nonslack_nodes,1)']),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_gnode-supply-flows.csv'],double([1:GN0;out.ss.gssol_all(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_gnode-demand-flows.csv'],double([1:GN0;out.ss.gdsol_all(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_slack-flows.csv'],double([out.ss.pn';out.ss.supp_flow(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-ratios.csv'],double([[1:out.ss.cn];out.ss.ccopt(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-discharge-pressure.csv'],double([[1:out.ss.cn];out.ss.csetopt(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-power.csv'],double([[1:out.ss.cn];out.ss.cpowopt(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_lmp-nodal-all.csv'],double([out.ss.fn';out.ss.trlmpnodal(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_lmp-gnodes.csv'],double([[1:GN0];out.ss.gnodelmp(1,:)]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ss_lmp-bidders.csv'],double([out.ss.gunique';out.ss.dglmp]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_pipe-lmp-in.csv'],double([pipe_cols;out.ss.lmpin(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_pipe-lmp-out.csv'],double([pipe_cols;out.ss.lmpout(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-lmp-in.csv'],double([1:out.ss.n0.nc;out.ss.lmpin(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-lmp-out.csv'],double([1:out.ss.n0.nc;out.ss.lmpout(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-pmax-mp.csv'],double([[1:out.ss.n0.nc];out.ss.mult0_pmax(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_comp-hpmax-mp.csv'],double([[1:out.ss.n0.nc];out.ss.mult0_cmax(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_flowbalrel.csv'],double([[1:out.ss.n0.nv];out.ss.flowbalrel(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_pipe-mass.csv'],double([pipe_cols;out.ss.pipe_mass_0(1,pipe_cols)]),'precision',16,'delimiter',',');
end
    if(par.out.intervals_out>0)
        %inputs on intervals
        out.ss.int_qbar=ss.int_qbar; out.ss.int_dmax=ss.int_dmax; out.ss.int_dmin=ss.int_dmin;
        out.ss.int_smax=ss.int_smax; out.ss.int_smin=ss.int_smin; out.ss.int_cd=ss.int_cd;
        out.ss.int_cs=ss.int_cs; out.ss.int_cslack=ss.int_cslack; out.ss.int_pslack=ss.int_pslack;
        if(ss.intervals>0 && ss.units==1)
            out.ss.int_cd=out.ss.int_cd*mmscfd_to_kgps;
            out.ss.int_cs=out.ss.int_cs*mmscfd_to_kgps; 
            out.ss.int_cslack=out.ss.int_cslack*mmscfd_to_kgps;
        end
        %------------------
        %revert to full gnode index list
        out.ss.gdsol_all=zeros(2,ss.n0.ng);
        out.ss.gdsol_all(:,ss.dmax_pos)=out.ss.gdsol;
        out.ss.gssol_all=zeros(2,ss.n0.ng);
        out.ss.gssol_all(:,ss.smax_pos)=out.ss.gssol;
        
        In=par.out.intervals_out;  %24 intervals on optimization period (1 day)
        T=ss.c.T/3600;             %time in hours
        int_bounds=[0:T/In:T]; out.ss.int_bounds=int_bounds;
        [out.ss.dbase_int]=pts_to_int(ss.m.xd'/3600,out.ss.dbase,int_bounds');
        [out.ss.gsub_int]=pts_to_int(ss.m.xd'/3600,out.ss.gsub,int_bounds');
        [out.ss.gslb_int]=pts_to_int(ss.m.xd'/3600,out.ss.gslb,int_bounds');
        [out.ss.gdub_int]=pts_to_int(ss.m.xd'/3600,out.ss.gdub,int_bounds');
        [out.ss.gdlb_int]=pts_to_int(ss.m.xd'/3600,out.ss.gdlb,int_bounds');
        [out.ss.gdsol_int]=pts_to_int(out.ss.tt0,out.ss.gdsol_all,int_bounds');
        [out.ss.gssol_int]=pts_to_int(out.ss.tt0,out.ss.gssol_all,int_bounds');
        [out.ss.dgflows_int]=pts_to_int(out.ss.tt0,out.ss.dgflows_all,int_bounds');
        [out.ss.supp_flow_int]=pts_to_int(out.ss.tt0,out.ss.supp_flow,int_bounds');
        [out.ss.flows_all_int]=pts_to_int(out.ss.tt0,out.ss.flows_all',int_bounds');
        [out.ss.Prslack_int]=pts_to_int(ss.m.xd'/3600,ss.m.Prslack',int_bounds');
        [out.ss.Prs_int]=pts_to_int(ss.m.xd'/3600,ss.m.Prs',int_bounds');
        [out.ss.Prd_int]=pts_to_int(ss.m.xd'/3600,ss.m.Prd',int_bounds');
        %------------------
        [out.ss.ppinopt_int]=pts_to_int(out.ss.tt0,out.ss.ppinopt,int_bounds');
        [out.ss.ppoutopt_int]=pts_to_int(out.ss.tt0,out.ss.ppoutopt,int_bounds');
        [out.ss.qqinopt_int]=pts_to_int(out.ss.tt0,out.ss.qqinopt,int_bounds');
        [out.ss.qqoutopt_int]=pts_to_int(out.ss.tt0,out.ss.qqoutopt,int_bounds');
        [out.ss.ppoptnodal_int]=pts_to_int(out.ss.tt0,out.ss.ppoptnodal,int_bounds');
        [out.ss.dgflows_int]=pts_to_int(out.ss.tt0,out.ss.dgflows,int_bounds');
        [out.ss.supp_flow_int]=pts_to_int(out.ss.tt0,out.ss.supp_flow,int_bounds');
        [out.ss.ccopt_int]=pts_to_int(out.ss.tt0,out.ss.ccopt,int_bounds');
        [out.ss.csetopt_int]=pts_to_int(out.ss.tt0,out.ss.csetopt,int_bounds');
        [out.ss.cpowopt_int]=pts_to_int(out.ss.tt0,out.ss.cpowopt,int_bounds');
        [out.ss.trlmpnodal_int]=pts_to_int(out.ss.tt0,out.ss.trlmpnodal,int_bounds');
        [out.ss.dglmp_int]=pts_to_int(out.ss.tt0,out.ss.dglmp,int_bounds');
        [out.ss.lmpin_int]=pts_to_int(out.ss.tt0,out.ss.lmpin,int_bounds');
        [out.ss.lmpout_int]=pts_to_int(out.ss.tt0,out.ss.lmpout,int_bounds');
        [out.ss.mult0_pmax_int]=pts_to_int(out.ss.tt0,out.ss.mult0_pmax,int_bounds');
        [out.ss.mult0_cmax_int]=pts_to_int(out.ss.tt0,out.ss.mult0_cmax,int_bounds');
        [out.ss.flowbalrel_int]=pts_to_int(out.ss.tt0,out.ss.flowbalrel,int_bounds');
        [out.ss.pipe_mass_int]=pts_to_int(out.ss.tt0,out.ss.pipe_mass_0,int_bounds');
        pipe_cols=[1:out.ss.n0.ne-out.ss.n0.nc]; comp_cols=[out.ss.n0.ne-out.ss.n0.nc+1:out.ss.n0.ne];
        

        if(par.out.steadystateonly==0)
            %write files
            dlmwrite([mfolder '\output_int_pipe-pressure-in.csv'],double([pipe_cols;out.ss.ppinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-pressure-out.ss.csv'],double([pipe_cols;out.ss.ppoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-pressure-in.csv'],double([1:out.ss.n0.nc;out.ss.ppinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-pressure-out.ss.csv'],double([1:out.ss.n0.nc;out.ss.ppoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-flow-in.csv'],double([pipe_cols;out.ss.qqinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-flow-out.ss.csv'],double([pipe_cols;out.ss.qqoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-flow-in.csv'],double([1:out.ss.n0.nc;out.ss.qqinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-flow-out.ss.csv'],double([1:out.ss.n0.nc;out.ss.qqoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_nodal-pressure.csv'],double([[1:out.ss.n0.nv];out.ss.ppoptnodal_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-physical-withdrawals.csv'],double([ss.m.guniqueind';out.ss.dgflows_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-supply-flows.csv'],double([ss.n0.phys_node';out.ss.gssol_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-demand-flows.csv'],double([ss.n0.phys_node';out.ss.gdsol_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_slack-flows.csv'],double([out.ss.pn';out.ss.supp_flow_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-ratios.csv'],double([[1:out.ss.cn];out.ss.ccopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-discharge-pressure.csv'],double([[1:out.ss.cn];out.ss.csetopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-power.csv'],double([[1:out.ss.cn];out.ss.cpowopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_lmp-nodal-all.csv'],double([out.ss.fn';out.ss.trlmpnodal_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_lmp-bidders.csv'],double([out.ss.gunique';out.ss.dglmp_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-lmp-in.csv'],double([pipe_cols;out.ss.lmpin_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-lmp-out.ss.csv'],double([pipe_cols;out.ss.lmpout_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-lmp-in.csv'],double([1:out.ss.n0.nc;out.ss.lmpin_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-lmp-out.ss.csv'],double([1:out.ss.n0.nc;out.ss.lmpout_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-pmax-mp.csv'],double([[1:out.ss.n0.nc];out.ss.mult0_pmax_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-hpmax-mp.csv'],double([[1:out.ss.n0.nc];out.ss.mult0_cmax_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_flowbalrel.csv'],double([[1:out.ss.n0.nv];out.ss.flowbalrel_int]),'precision',16,'delimiter',',');
        end
    end

if(ss.m.save_state==1)
dlmwrite([mfolder '\output_ss_state_save.csv'],double(full(ss.state_save)),'precision',16,'delimiter',',');
end

par.out.ss=out.ss;
%par.ss=ss;

%out.ss.mult0_pmax=ss.mult0_pmax/2*ss.m.N/(ss.c.psc/1000000)/mpa_to_psi;    %output pressure marginal prices ($/
%out.ss.mult0_cmax=ss.mult0_cmax/2*ss.m.N*3.6/0.75; %compression marginal prices ($/hp)

function [xints]=pts_to_int(tpts,xpts,ibnds)
    xbnds=interp1qr(tpts,xpts,ibnds); In=length(ibnds)-1;
    xints=(xbnds(1:In,:)+xbnds(2:In+1,:))/2;
return;