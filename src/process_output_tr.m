function [par]=process_output_tr(par)
% Anatoly Zlotnik, January 2019

mfolder=par.mfolder;
out=par.out;

if(par.out.dosim==1), sim=par.sim; end
tr=par.tr; 
psi_to_pascal=tr.c.psi_to_pascal;
mpa_to_psi=1000000/psi_to_pascal;
tr.c.mpa_to_psi=mpa_to_psi;
mmscfd_to_kgps=tr.c.mmscfd_to_kgps;
if(par.out.doZ==1), b1=tr.c.b1; b2=tr.c.b2; end

%process simulation output
if(par.out.dosim==1)
out.tt=sim.tt/3600;         %plotting time in hours
out.qqinsim=sim.qqin;                  %flow boundary in
out.qqoutsim=sim.qqout;                %flow boundary out
if(par.out.plotnodal==1)
    out.ppsim=sim.pnodout(:,1:sim.n0.nv);  %pressure (nodal)
    out.qqsim=[out.qqinsim out.qqoutsim];  %all boundary flows
else
    out.ppsim=sim.pnodout(:,1:sim.n.nv);   %pressure (all)
    out.qqsim=sim.qq;                      %flows (all)
end
if(par.out.units==1), out.ppsim=out.ppsim/psi_to_pascal;
out.qqsim=out.qqsim/mmscfd_to_kgps;
out.qqinsim=out.qqinsim/mmscfd_to_kgps;  out.qqoutsim=out.qqoutsim/mmscfd_to_kgps;end
end

%process optimization output
%process dimensional solution (metric)
pp=zeros(length(tr.tt0),tr.n0.nv);
%s_nodal=[tr.m.pslack tr.m.pslack(:,1)];
s_nodal=tr.m.pslack;
if(par.tr.m.extension>0), s_nodal=tr.m.pslack; end
scf=find(ismember(tr.m.snodes,tr.m.comp_pos(tr.m.spos,1)));
s(scf,:)=tr.cc0(tr.m.spos,:).*s_nodal(scf,:);
pp(:,tr.m.snodes)=s';
pp(:,tr.m.dnodes)=tr.pp0';
qq=tr.qq0'*tr.m.Xs*tr.c.qsc; 
ff=tr.qq0'*tr.c.qsc; 
for j=1:length(tr.m.dpos)
    cpj=tr.n.comp_pos(tr.m.dpos(j),1);
    pp(:,cpj)=pp(:,cpj).*tr.cc0(tr.m.dpos(j),:)';
end
%nodal pressure (before compressors)
pp_nodal=zeros(length(tr.tt0),tr.n.nv);
pp_nodal(:,tr.m.snodes)=s_nodal';
pp_nodal(:,tr.m.dnodes)=tr.pp0';

if(tr.m.doZ==1), pp=rho_to_p_nd(pp,tr.c.b1,tr.c.b2,tr.c.psc);
    pp_nodal=rho_to_p_nd(pp_nodal,tr.c.b1,tr.c.b2,tr.c.psc); end

%compute mass in pipes
if(tr.m.doZ==1), p_density=p_to_rho(pp',tr.c.b1,tr.c.b2,par.tr);
else p_density=tr.c.psc*pp'/(tr.c.gasR*tr.c.gasT); end
p_comp=par.tr.cc0;
out.p_mass=pipe_mass(p_density,p_comp,par.tr.m);    %all pipes
out.pipe_mass_0=(par.tr.n.disc_to_edge*out.p_mass)';       %original pipes

tr.pnodin=pp_nodal*tr.c.psc;   %nodal pressures (before compression)
tr.pnodout=pp*tr.c.psc;     %nodal pressures (after compression)
tr.qqin=qq(:,tr.n.from_flows);
tr.qqout=qq(:,tr.n.to_flows);
%pipe inlet and outlet pressure (compressors only at inlets)
tr.ppin=zeros(par.tr.m.N1+1,tr.n0.ne);
tr.ppout=zeros(par.tr.m.N1+1,tr.n0.ne);
for j=1:tr.n0.ne
    if(tr.n0.comp_bool(j)==1)
        %tr.ppin(:,j)=tr.pnodout(:,tr.n0.from_id(j));
        tr.ppin(:,j)=tr.pnodin(:,tr.n0.from_id(j));
        tr.ppout(:,j)=tr.pnodin(:,tr.n0.to_id(j));
    elseif(tr.n0.comp_bool(j)==0)
        tr.ppin(:,j)=tr.pnodin(:,tr.n0.from_id(j));
        tr.ppout(:,j)=tr.pnodin(:,tr.n0.to_id(j));
    end
end
%pipe inlet and outlet flux
tr.ffin=ff(:,tr.n.from_flows);
tr.ffout=ff(:,tr.n.to_flows);

out.tt0=tr.tt0/3600;         %plotting time in hours
out.qqinopt=tr.qqin;                  %flow boundary in
out.qqoutopt=tr.qqout;                %flow boundary out
%if(par.out.plotnodal==1)
    out.ppoptnodal=tr.pnodin(:,1:tr.n0.nv);
    out.ppopt=tr.pnodout(:,1:tr.n0.nv);  %pressure (nodal)
    out.qqopt=[out.qqinopt out.qqoutopt];  %all boundary flows
%else
    %out.ppoptall=tr.pnodout(:,1:tr.n.nv);   %pressure (all)
    %out.qqoptall=tr.qq0;                      %flows (all)
%end
    out.ppinopt=tr.ppin; out.ppoutopt=tr.ppout;
if(par.out.units==1), out.ppopt=out.ppopt/psi_to_pascal; out.ppoptnodal=out.ppoptnodal/psi_to_pascal; 
out.ppinopt=tr.ppin/psi_to_pascal; out.ppoutopt=tr.ppout/psi_to_pascal; 
out.qqopt=out.qqopt/mmscfd_to_kgps; out.qqinopt=out.qqinopt/mmscfd_to_kgps;  
out.qqoutopt=out.qqoutopt/mmscfd_to_kgps; out.pipe_mass_0=out.pipe_mass_0/mmscfd_to_kgps/86400; end

%market flow solution
tr.m.Yd=interp1qr(tr.m.xd',tr.m.Yq(1:tr.m.FN,:)',tr.tt0)';
tr.m.Yd(tr.m.guniqueind,:)=tr.m.Yd(tr.m.guniqueind,:)+tr.m.gtod*tr.fd0;
tr.m.Ygd=tr.fd0(1:length(tr.m.gd),:);
tr.m.Ygs=-tr.fd0(length(tr.m.gd)+1:length(tr.m.gall),:);

%process comparison
if(par.out.dosim==1)
out.ttcom=min((sim.periodgrid'*tr.c.Tsc)/3600,out.tt0(end)-eps);
ttcomsim=min((sim.periodgrid'*tr.c.Tsc+(sim.nperiods-1)*tr.c.T)/3600,max(out.tt)-eps);
psim1=interp1qr(out.tt,out.ppsim,ttcomsim);popt1=interp1qr(out.tt0,out.ppopt,out.ttcom);
qsim1=interp1qr(out.tt,out.qqsim,ttcomsim);qopt1=interp1qr(out.tt0,out.qqopt,out.ttcom);
out.pdiff=psim1-popt1; out.qdiff=qsim1-qopt1; 
%out.prel=out.pdiff./psim1; out.qrel=out.qdiff/mean(mean(abs(sim.qq)));
out.prel=2*out.pdiff./(psim1+popt1); out.qrel=2*out.qdiff/mean(mean(abs(qsim1)+abs(qopt1)));
%out.prel=out.pdiff./(psim1+popt1)*2; out.qrel=out.qdiff./(qsim1+qopt1)*2;
end

%compressor discharge pressures
if(par.out.dosim==1), out.csetsim=psim1(:,tr.m.comp_pos(:,1)); end
out.csetopt=out.ppopt(:,tr.m.comp_pos(:,1));

%process parameters (compression ratios and demands)
out.cc=tr.cc0'; out.td=tr.m.xd/3600; 
out.dbase=tr.m.Yq(1:length(tr.m.fn),:)'*tr.c.qsc;     %base flow "q"
out.gsub=tr.m.Yubs'*tr.c.qsc;   %upper bounds on sales
out.gslb=tr.m.Ylbs'*tr.c.qsc;   %lower bounds on sales
out.gdub=tr.m.Yubd'*tr.c.qsc;   %upper bounds on buys
out.gdlb=tr.m.Ylbd'*tr.c.qsc;   %lower bounds on buys
out.gdsol=tr.m.Ygd'*tr.c.qsc;   %gnode buyer solutions
out.gssol=-tr.m.Ygs'*tr.c.qsc;   %gnode seller solutions
%gnode buyer and seller solutions for all original gnodes 
GN0=length(tr.n0.phys_node);    %number of original gnodes
out.gdsol_all=zeros(length(tr.tt0),GN0); out.gdsol_all(:,tr.dmax_pos)=out.gdsol;
out.gssol_all=zeros(length(tr.tt0),GN0); out.gssol_all(:,tr.smax_pos)=out.gssol;
out.dgflows=full(tr.m.Yd(tr.m.guniqueind,:))'*tr.c.qsc; %flow at nodes with gnodes
%out.dgflows_all=full(tr.m.Yd(tr.m.guniqueind,:))'*tr.c.qsc; %flow at all original gnodes
slinks=tr.m.comp_pos(tr.m.spos,2); 
out.supp_flow=qq(:,slinks);    %supply flow 
out.nonslack_flow=full(tr.m.Yd(1:length(tr.n0.nonslack_nodes),:))*tr.c.qsc;
out.flows_all=zeros(tr.n0.nv,par.tr.m.N1+1);
out.flows_all(tr.n0.slack_nodes,:)=-out.supp_flow;
out.flows_all(tr.n0.nonslack_nodes,:)=out.nonslack_flow;
out.dgflows_all=out.flows_all(tr.n0.phys_node,:)'; %flow at all original gnodes
% out.supp_flow_sim=sim.qq(:,slinks);    %supply flow 
% out.flows_all_sim=zeros(sim.n0.nv,length(sim.tt));
% out.flows_all_sim(sim.n0.slack_nodes,:)=-out.supp_flow_sim;
% out.flows_all_sim(sim.n0.nonslack_nodes,:)=full(sim.m.Yd(1:length(sim.n0.nonslack_nodes),:))*sim.c.qsc;
if(par.out.units==1), out.dbase=out.dbase/mmscfd_to_kgps; out.gsub=out.gsub/mmscfd_to_kgps; 
    out.gslb=out.gslb/mmscfd_to_kgps; out.gdub=out.gdub/mmscfd_to_kgps; out.gdlb=out.gdlb/mmscfd_to_kgps;
    out.gdsol=out.gdsol/mmscfd_to_kgps; out.gssol=out.gssol/mmscfd_to_kgps;
    out.gdsol_all=out.gdsol_all/mmscfd_to_kgps; out.gssol_all=out.gssol_all/mmscfd_to_kgps;
    out.dgflows=out.dgflows/mmscfd_to_kgps; out.dgflows_all=out.dgflows_all/mmscfd_to_kgps; 
    out.supp_flow=out.supp_flow/mmscfd_to_kgps; 
    out.flows_all=out.flows_all/mmscfd_to_kgps;  
    %out.supp_flow_sim/mmscfd_to_kgps; %out.flows_all_sim/mmscfd_to_kgps;  
end


%check nodal flow balance
out.flowbal=tr.n0.Amp*out.qqoutopt'+tr.n0.Amm*out.qqinopt'-out.flows_all;
out.flowbalrel=3*out.flowbal./(abs(tr.n0.Amp*out.qqoutopt')+abs(tr.n0.Amm*out.qqinopt')+abs(out.flows_all));
out.flowbalrel(mean(out.flowbal')./mean(out.flowbalrel')<tr.m.opt_tol,:)=0; out.flowbalrel=out.flowbalrel';

%out.flowbals=tr.n0.Amp*out.qqoutsim'+tr.n0.Amm*out.qqinsim'-out.flows_all;
%out.flowbalsrel=3*out.flowbal./(abs(tr.n0.Amp*out.qqoutsim')+abs(tr.n0.Amm*out.qqinsim')+abs(out.flows_all));
%out.flowbalsrel(mean(out.flowbal')./mean(out.flowbalrel')<tr.m.opt_tol,:)=0;

%compressor power
if(par.out.dosim==1)
    cpossim=par.tr.m.comp_pos; m=tr.m.mpow;
    out.cccom=interp1qr(out.tt0,tr.cc0',out.ttcom);
    qcompsim=interp1qr(out.tt,sim.qq(:,cpossim(:,2)),ttcomsim); 
    cpow_nd=(abs(qcompsim)).*((out.cccom).^(m)-1);
    out.cpowsim=cpow_nd.*kron(tr.m.eff',ones(size(cpow_nd,1),1))*tr.c.mmscfd_to_hp/mmscfd_to_kgps;
end
out.ccopt=tr.cc0';
cposopt=par.tr.m.comp_pos; m=tr.m.mpow;
qcompopt=qq(:,cposopt(:,2)); cpow_nd=(abs(qcompopt)).*((tr.cc0').^(m)-1);
out.cpowopt=cpow_nd.*kron(tr.m.eff',ones(size(cpow_nd,1),1))*tr.c.mmscfd_to_hp/mmscfd_to_kgps;

%process locational marginal price
% out.lmptr=par.tr.lmp0(par.tr.m.flexnodes,:)'/2*par.tr.m.N/par.tr.m.odw*par.tr.c.Tsc*par.tr.c.Tsc/2;
% lmpss=par.tr.lmp0(par.tr.m.flexnodes,:)'/par.tr.m.odw*par.tr.c.Tsc*par.tr.c.Tsc/2;
% out.lmptr=par.tr.lmp0'/2*par.tr.m.N/par.tr.m.odw*par.tr.c.Tsc*par.tr.c.Tsc/2;
% lmpss=par.tr.lmp0'/tr.m.odw*tr.c.Tsc*tr.c.Tsc/2;
if(tr.m.N==0), trmN=2; end
if(tr.m.N>1), trmN=tr.m.N; end
out.trlmp=tr.lmp0'/2*trmN;    %all lmps
out.trlmpnodal=out.trlmp(:,1:length(tr.m.fn));
out.trlmpnodal_all=zeros(tr.m.N1+1,tr.n0.nv);
out.trlmpnodal_all(:,tr.n0.slack_nodes)=-tr.m.prslack';
out.trlmpnodal_all(:,tr.n0.nonslack_nodes)=out.trlmpnodal;
out.gnodelmp=out.trlmpnodal_all(:,tr.n0.phys_node);
%out.gnodelmp=tr.lmp0(tr.n0.phys_node,:)'/2*trmN;     %lmps at all gnodes
out.dglmp=tr.lmp0(tr.m.guniqueind,:)'/2*trmN;     %lmps at market nodes
out.gdlmp=tr.lmp0(tr.m.gallind(1:length(tr.m.gd)),:)'/2*trmN;     %lmps at demand gnodes
out.gslmp=tr.lmp0(tr.m.gallind(length(tr.m.gd):length(tr.m.gall)),:)'/2*trmN;     %lmps at supply gnodes
% if(par.out.dosim==1), out.lmpss=par.tr.lmp0'; end
out.Prd=tr.m.Prd; out.Prs=tr.m.Prs; 
out.Prslack=interp1qr(tr.m.xd',tr.m.Prslack',out.tt0);  %bid and offer prices
out.mult0_pmax=tr.mult0_pmax'/2*trmN*tr.c.psi_to_pascal/tr.c.psc*3600;    %output pressure marginal prices ($/Psi/hr)
out.mult0_cmax=tr.mult0_cmax'/2*trmN*3.6/0.75; %compression marginal prices ($/hp)
if(par.out.units==1), out.trlmp=out.trlmp*mmscfd_to_kgps; 
     out.trlmpnodal=out.trlmpnodal*mmscfd_to_kgps; 
     out.dglmp=out.dglmp*mmscfd_to_kgps; out.gnodelmp=out.gnodelmp*mmscfd_to_kgps;
     out.gdlmp=out.gdlmp*mmscfd_to_kgps; out.gslmp=out.gslmp*mmscfd_to_kgps; 
     out.Prd=out.Prd*mmscfd_to_kgps; out.Prs=out.Prs*mmscfd_to_kgps;
     out.Prslack=out.Prslack*mmscfd_to_kgps;
     %if(par.out.dosim==1), out.lmpss=out.lmpss*mmscfd_to_kgps; end
end
out.lmpin=zeros(size(out.ppinopt)); out.lmpout=zeros(size(out.ppoutopt));
for j=1:tr.n0.ne
    if(ismember(tr.n0.from_id(j),tr.m.pn))
        out.lmpin(:,j)=out.Prslack(:,find(tr.m.pn==tr.n0.from_id(j))); else
        out.lmpin(:,j)=out.trlmpnodal(:,find(tr.m.fn==tr.n0.from_id(j))); end
    if(ismember(tr.n0.to_id(j),tr.m.pn))
        out.lmpout(:,j)=out.Prslack(:,find(tr.m.pn==tr.n0.to_id(j))); else
        out.lmpout(:,j)=out.trlmpnodal(:,find(tr.m.fn==tr.n0.to_id(j))); end
end
%cmap=colormap;
%set(0,'DefaultAxesColorOrder',cmap(floor(rand(length(tr.m.gs),1)*64)+1,:))

out.guniqueind=tr.m.guniqueind; out.gunique=tr.m.gunique; out.fn=tr.m.fn; out.pn=tr.m.pn;
out.n0=tr.n0; out.n=tr.n; out.gd=tr.m.gd; out.gs=tr.m.gs; out.FN=tr.m.FN; out.PN=tr.m.PN;
out.cn=tr.m.C; 
out.mfolder=par.mfolder;

if(par.out.savecsvoutput==1)
    if(par.out.intervals_out==0)
        pipe_cols=[1:out.n0.ne-out.n0.nc]; comp_cols=[out.n0.ne-out.n0.nc+1:out.n0.ne];
        dlmwrite([mfolder '\output_ts_tpts.csv'],double(out.tt0),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_pipe-pressure-in.csv'],double([pipe_cols;out.ppinopt(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_pipe-pressure-out.csv'],double([pipe_cols;out.ppoutopt(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-pressure-in.csv'],double([1:out.n0.nc;out.ppinopt(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-pressure-out.csv'],double([1:out.n0.nc;out.ppoutopt(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_pipe-flow-in.csv'],double([pipe_cols;out.qqinopt(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_pipe-flow-out.csv'],double([pipe_cols;out.qqoutopt(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-flow-in.csv'],double([1:out.n0.nc;out.qqinopt(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-flow-out.csv'],double([1:out.n0.nc;out.qqoutopt(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_nodal-pressure.csv'],double([[1:out.n0.nv];out.ppoptnodal]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ts_gnode-physical-withdrawals.csv'],double([out.gunique';out.dgflows]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_nonslack-flows.csv'],double([out.fn';out.nonslack_flow']),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_gnode-supply-flows.csv'],double([[1:GN0];out.gssol_all]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_gnode-demand-flows.csv'],double([[1:GN0];out.gdsol_all]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_slack-flows.csv'],double([out.pn';out.supp_flow]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-ratios.csv'],double([[1:out.cn];out.ccopt]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-discharge-pressure.csv'],double([[1:out.cn];out.csetopt]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-power.csv'],double([[1:out.cn];out.cpowopt]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_lmp-nodal-all.csv'],double([out.fn';out.trlmpnodal]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_lmp-gnodes.csv'],double([[1:GN0];out.gnodelmp]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ts_lmp-bidders.csv'],double([out.gunique';out.dglmp]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_pipe-lmp-in.csv'],double([pipe_cols;out.lmpin(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_pipe-lmp-out.csv'],double([pipe_cols;out.lmpout(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-lmp-in.csv'],double([1:out.n0.nc;out.lmpin(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-lmp-out.csv'],double([1:out.n0.nc;out.lmpout(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-pmax-mp.csv'],double([[1:out.n0.nc];out.mult0_pmax]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_comp-hpmax-mp.csv'],double([[1:out.n0.nc];out.mult0_cmax]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_flowbalrel.csv'],double([[1:out.n0.nv];out.flowbalrel]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ts_pipe-mass.csv'],double([pipe_cols;out.pipe_mass_0(:,pipe_cols)]),'precision',16,'delimiter',',');
    end
    if(par.out.intervals_out>0)
        %inputs on intervals
        out.int_qbar=tr.int_qbar; out.int_dmax=tr.int_dmax; out.int_dmin=tr.int_dmin;
        out.int_smax=tr.int_smax; out.int_smin=tr.int_smin; out.int_cd=tr.int_cd;
        out.int_cs=tr.int_cs; out.int_cslack=tr.int_cslack; out.int_pslack=tr.int_pslack;
        if(tr.intervals>0 && tr.units==1)
            out.int_cd=out.int_cd*mmscfd_to_kgps;
            out.int_cs=out.int_cs*mmscfd_to_kgps; 
            out.int_cslack=out.int_cslack*mmscfd_to_kgps;
        end
        %------------------
        %revert to full gnode index list
        out.gdsol_all=zeros(par.out.intervals_out+1,tr.n0.ng);
        out.gdsol_all(:,tr.dmax_pos)=out.gdsol;
        out.gssol_all=zeros(par.out.intervals_out+1,tr.n0.ng);
        out.gssol_all(:,tr.smax_pos)=out.gssol;
        
        In=par.out.intervals_out;  %24 intervals on optimization period (1 day)
        T=tr.c.T/3600;             %time in hours
        int_bounds=[0:T/In:T]; out.int_bounds=int_bounds;
        [out.dbase_int]=pts_to_int(tr.m.xd'/3600,out.dbase,int_bounds');
        [out.gsub_int]=pts_to_int(tr.m.xd'/3600,out.gsub,int_bounds');
        [out.gslb_int]=pts_to_int(tr.m.xd'/3600,out.gslb,int_bounds');
        [out.gdub_int]=pts_to_int(tr.m.xd'/3600,out.gdub,int_bounds');
        [out.gdlb_int]=pts_to_int(tr.m.xd'/3600,out.gdlb,int_bounds');
        [out.gdsol_int]=pts_to_int(out.tt0,out.gdsol_all,int_bounds');
        [out.gssol_int]=pts_to_int(out.tt0,out.gssol_all,int_bounds');
        [out.dgflows_int]=pts_to_int(out.tt0,out.dgflows_all,int_bounds');
        [out.supp_flow_int]=pts_to_int(out.tt0,out.supp_flow,int_bounds');
        [out.nonslack_flow_int]=pts_to_int(out.tt0,out.nonslack_flow',int_bounds');
        [out.flows_all_int]=pts_to_int(out.tt0,out.flows_all',int_bounds');
        [out.Prslack_int]=pts_to_int(tr.m.xd'/3600,tr.m.Prslack',int_bounds');
        [out.Prs_int]=pts_to_int(tr.m.xd'/3600,tr.m.Prs',int_bounds');
        [out.Prd_int]=pts_to_int(tr.m.xd'/3600,tr.m.Prd',int_bounds');
        %------------------
        [out.ppinopt_int]=pts_to_int(out.tt0,out.ppinopt,int_bounds');
        [out.ppoutopt_int]=pts_to_int(out.tt0,out.ppoutopt,int_bounds');
        [out.qqinopt_int]=pts_to_int(out.tt0,out.qqinopt,int_bounds');
        [out.qqoutopt_int]=pts_to_int(out.tt0,out.qqoutopt,int_bounds');
        [out.ppoptnodal_int]=pts_to_int(out.tt0,out.ppoptnodal,int_bounds');
        [out.dgflows_int]=pts_to_int(out.tt0,out.dgflows,int_bounds');
        [out.supp_flow_int]=pts_to_int(out.tt0,out.supp_flow,int_bounds');
        [out.ccopt_int]=pts_to_int(out.tt0,out.ccopt,int_bounds');
        [out.csetopt_int]=pts_to_int(out.tt0,out.csetopt,int_bounds');
        [out.cpowopt_int]=pts_to_int(out.tt0,out.cpowopt,int_bounds');
        [out.trlmpnodal_int]=pts_to_int(out.tt0,out.trlmpnodal,int_bounds');
        [out.gnodelmp_int]=pts_to_int(out.tt0,out.gnodelmp,int_bounds');
        [out.dglmp_int]=pts_to_int(out.tt0,out.dglmp,int_bounds');
        [out.lmpin_int]=pts_to_int(out.tt0,out.lmpin,int_bounds');
        [out.lmpout_int]=pts_to_int(out.tt0,out.lmpout,int_bounds');
        [out.mult0_pmax_int]=pts_to_int(out.tt0,out.mult0_pmax,int_bounds');
        [out.mult0_cmax_int]=pts_to_int(out.tt0,out.mult0_cmax,int_bounds');
        [out.flowbalrel_int]=pts_to_int(out.tt0,out.flowbalrel,int_bounds');
        [out.pipe_mass_int]=pts_to_int(out.tt0,out.pipe_mass_0,int_bounds');
        pipe_cols=[1:out.n0.ne-out.n0.nc]; comp_cols=[out.n0.ne-out.n0.nc+1:out.n0.ne];
        
        %write files
        dlmwrite([mfolder '\output_int_pipe-pressure-in.csv'],double([pipe_cols;out.ppinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_pipe-pressure-out.csv'],double([pipe_cols;out.ppoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-pressure-in.csv'],double([1:out.n0.nc;out.ppinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-pressure-out.csv'],double([1:out.n0.nc;out.ppoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_pipe-flow-in.csv'],double([pipe_cols;out.qqinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_pipe-flow-out.csv'],double([pipe_cols;out.qqoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-flow-in.csv'],double([1:out.n0.nc;out.qqinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-flow-out.csv'],double([1:out.n0.nc;out.qqoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_nodal-pressure.csv'],double([[1:out.n0.nv];out.ppoptnodal_int]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_int_gnode-physical-withdrawals.csv'],double([tr.m.guniqueind';out.dgflows_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_nonslack-flows.csv'],double([out.fn';out.nonslack_flow_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_gnode-supply-flows.csv'],double([1:GN0;out.gssol_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_gnode-demand-flows.csv'],double([1:GN0;out.gdsol_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_slack-flows.csv'],double([out.pn';out.supp_flow_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-ratios.csv'],double([[1:out.cn];out.ccopt_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-discharge-pressure.csv'],double([[1:out.cn];out.csetopt_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-power.csv'],double([[1:out.cn];out.cpowopt_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_lmp-nodal-all.csv'],double([out.fn';out.trlmpnodal_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_lmp-gnodes.csv'],double([[1:GN0];out.gnodelmp_int]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_int_lmp-bidders.csv'],double([out.gunique';out.dglmp_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_pipe-lmp-in.csv'],double([pipe_cols;out.lmpin_int(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_pipe-lmp-out.csv'],double([pipe_cols;out.lmpout_int(:,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-lmp-in.csv'],double([1:out.n0.nc;out.lmpin_int(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-lmp-out.csv'],double([1:out.n0.nc;out.lmpout_int(:,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-pmax-mp.csv'],double([[1:out.n0.nc];out.mult0_pmax_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_comp-hpmax-mp.csv'],double([[1:out.n0.nc];out.mult0_cmax_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_flowbalrel.csv'],double([[1:out.n0.nv];out.flowbalrel_int]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_int_pipe-mass.csv'],double([pipe_cols;out.pipe_mass_int(:,pipe_cols)]),'precision',16,'delimiter',',');
    end
end
if(tr.m.save_state==1)
dlmwrite([mfolder '\output_tr_state_save.csv'],double(full(tr.state_save)),'precision',16,'delimiter',',');
end

par.out=out;
par.tr=tr;

%out.mult0_pmax=tr.mult0_pmax/2*tr.m.N/(tr.c.psc/1000000)/mpa_to_psi;    %output pressure marginal prices ($/
%out.mult0_cmax=tr.mult0_cmax/2*tr.m.N*3.6/0.75; %compression marginal prices ($/hp)

function [xints]=pts_to_int(tpts,xpts,ibnds)
    xbnds=interp1qr(tpts,xpts,ibnds); In=length(ibnds)-1;
    xints=(xbnds(1:In,:)+xbnds(2:In+1,:))/2;
return;