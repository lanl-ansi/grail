function []=gas_out_plots_i(par)

out=par.out; mfolder=out.ss.mfolder;

if(par.out.ss.plotlarge==1), 
    plotpos=[20,100,1600,860]; paperpos=[0 0 1600 860]; papersize=[1600 860]; else
    plotpos=[20,400,1200,400]; paperpos=[0 0 1200 400]; papersize=[1200 400]; end

int_bounds=out.ss.int_bounds;
tp=kron(int_bounds,[1 1])';tp(1)=[];tp(end)=[];

if(par.out.ss.plotmarketflowlims==1)
    %fixed gas withdrawals, bounds on flexible demands and supplies
    f1=figure(1); clf
    set(f1,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    tv=kron(out.ss.int_qbar(:,out.ss.guniqueind)',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    if(par.out.ss.units==1), title('Baseline Gas Withdrawals (mmscfd)','fontweight','bold'), else
        title('Baseline Gas Withdrawals (kg/s)','fontweight','bold'), end
    if(length(out.ss.guniqueind)<20), legend(num2str(out.ss.gunique),'Location','SouthEast'), end
    if(~isempty(par.ss.dmax_pos))
        subaxis(1,3,2,'SpacingHoriz',0.05), 
        tv=kron(out.ss.int_dmax',[1 1])';
        plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), hold on,
        tv=kron(out.ss.int_dmin',[1 1])';
        plot(tp,tv,'--','LineWidth',3), hold off, legend(num2str(par.ss.dmax_pos),'Location','SouthEast')
        if(par.out.ss.units==1), title('Bounds on Buyer Offtakes (mmscfd)','fontweight','bold'), else
        title('Bounds on Buyer Offtakes (kg/s))','FontWeight','bold'), end
    end
    if(~isempty(par.ss.smax_pos))
        subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05), 
        tv=kron(out.ss.int_smax',[1 1])';
        plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), hold on,
        tv=kron(out.ss.int_smin',[1 1])';
        plot(tp,tv,'--','LineWidth',3), hold off, legend(num2str(par.ss.smax_pos),'Location','SouthEast')
        if(par.out.ss.units==1), title('Bounds on Seller Supplies (mmscfd)','fontweight','bold'), else
        title('Bounds on Seller Supplies (kg/s))','FontWeight','bold'), end
    end
    if(par.out.ss.plotpdf==1)
        set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
            'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
        eval(['print -dpdf ' mfolder '\1marketflowlims.pdf']), end
    if(par.out.ss.ploteps==1)
        set(gcf,'PaperPositionMode','auto')
        eval(['print -depsc ' mfolder '\1marketflowlims.eps'])
    end
end

if(par.out.ss.plotmarketpricebids==1)
    %price bids and offers
    f2=figure(2); clf
    set(f2,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    tv=kron(out.ss.int_cslack',[1 1])';
	plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    if(par.out.ss.units==1), title('Pressure (slack) Node Prices ($/mscf)','fontweight','bold'), else
    title('Pressure (slack) Node (10$/kg)','FontWeight','bold'), end
    legend(num2str(out.ss.PN),'Location','SouthEast')
    if(~isempty(par.ss.dmax_pos))
        subaxis(1,3,2,'SpacingHoriz',0.05),
        tv=kron(out.ss.int_cd',[1 1])';
        plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
        if(par.out.ss.units==1), title('Demand Bid Prices ($/mscf)','fontweight','bold'), else
        title('Demand Bid Prices (10$/kg)','FontWeight','bold'), end
        legend(num2str(par.ss.dmax_pos),'Location','SouthEast')
    end
    if(~isempty(par.ss.smax_pos))
        subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05), 
        tv=kron(out.ss.int_cs',[1 1])';
        plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
        if(par.out.ss.units==1), title('Supply Offer Prices ($/mscf)','fontweight','bold'), else
        title('Supply Offer Prices (10$/kg)','FontWeight','bold'), end
        legend(num2str(par.ss.smax_pos),'Location','SouthEast')
    end
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\2marketpricebids.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\2marketpricebids.eps']), end
end


if(par.out.ss.plotsolflow==1)
    %flow schedule solution
    f3=figure(3); clf
    set(f3,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    tv=kron([-out.ss.supp_flow_int out.ss.dgflows_int]',[1 1])';
    plot(tp,tv,'LineWidth',3), xlabel('hours'), ylabel(''), hold off
    if(par.out.ss.units==1), title('Cleared Nodal Gas Withdrawals (mmscfd)','fontweight','bold'), else
        title('Cleared Nodal Gas Withdrawals (kg/s)','fontweight','bold'), end
    legend(num2str([out.ss.PN;out.ss.gunique]),'Location','SouthEast')
    subaxis(1,3,2,'SpacingHoriz',0.05),
    tv=kron(out.ss.gdsol_int',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    if(par.out.ss.units==1), title('Demand gNode Purchases (mmscfd)','fontweight','bold'), else
        title('Demand gNode Purchases  (kg/s)','fontweight','bold'), end
    legend(num2str(out.ss.gd),'Location','SouthEast')
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05),
    tv=kron([out.ss.supp_flow_int out.ss.gssol_int]',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    if(par.out.ss.units==1), title('Supply gNode Sales (mmscfd)','fontweight','bold'), else
        title('Supply gNode Sales  (kg/s)','fontweight','bold'), end
    legend(num2str([out.ss.PN;out.ss.gs]),'Location','SouthEast')
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\3solflow.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\3solflow.eps']), end
end

if(par.out.ss.plotsolprice==1)
    %flexible gas withdrawals, lmps
    f4=figure(4); clf
    set(f4,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,2,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    tv=kron(out.ss.trlmpnodal_int',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    if(par.out.ss.units==1), title('All Nodal LMPs ($/mscf)','fontweight','bold'), else
    title('All Nodal LMPs (10$/kg)','FontWeight','bold'), end
    subaxis(1,2,2,'MarginRight',0.05,'SpacingHoriz',0.05),
    tv=kron(out.ss.dglmp_int',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    if(par.out.ss.units==1), title('gNode LMPs ($/mscf)','fontweight','bold'), else
    title('gNode LMPs (10$/kg)','FontWeight','bold'), end
    legend(num2str(out.ss.gunique),'Location','SouthEast')
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\4solprice.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\4solprice.eps']), end
end

if(par.out.ss.plotopt==1)
    %optimization plot
    f5=figure(5);
    set(f5,'position',[20,400,800,400],'Color',[1 1 1]);
    subplot(2,1,1), plot(out.ss.tt0,out.ss.ppopt,'LineWidth',3), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Optimization Solution Nodal Pressures (psi)','fontweight','bold'), else
        title('Optimization Solution Nodal Pressures (MPa)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.ss.tt0,out.ss.qqopt,'LineWidth',3), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Optimization Solution Flows (mmscfd)','fontweight','bold'), else
        title('Optimization Solution Flows (kg/s)','fontweight','bold'), end   
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\5opt.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\5opt.eps']), end
end

if(par.out.ss.dosim==1)
if(par.out.ss.plotsim==1)
    %simulation plot
    f6=figure(6);
    set(f6,'position',[20,400,800,400],'Color',[1 1 1]);
    subplot(2,1,1), plot(out.ss.tt,out.ss.ppsim), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Simulation Solution Nodal Pressures (psi)','fontweight','bold'), else
        title('Simulation Solution Nodal Pressures (MPa)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.ss.tt,out.ss.qqsim), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Simulation Solution Flows (mmscfd)','fontweight','bold'), else
        title('Simulation Solution Flows (kg/s)','fontweight','bold'), end
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\6sim.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\6sim.eps']), end
end

if(par.out.ss.plotcompabs==1)
    %comparison plot absolute
    f7=figure(7);
    set(f7,'position',[20,400,800,400],'Color',[1 1 1]);
    if(par.out.ss.plotabsolute==1), pdp=abs(out.ss.pdiff); else, pdp=out.ss.pdiff; end
    if(par.out.ss.plotabsolute==1), pdq=abs(out.ss.qdiff); else, pdq=out.ss.qdiff; end
    subplot(2,1,1), plot(out.ss.ttcom,pdp), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Pressure Absolute Difference (psi)','fontweight','bold'), else
        title('Pressure Absolute Difference (MPa)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.ss.ttcom,pdq), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Flow Absolute Difference (mmscfd)','fontweight','bold'), else
        title('Flow Absolute Difference (kg/s)','fontweight','bold'), end   
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\7diff.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\7diff.eps']), end
end

if(par.out.ss.plotcomprel==1)
    %comparison plot relative
    f8=figure(8);
    set(f8,'position',[20,400,800,400],'Color',[1 1 1]);
    if(par.out.ss.plotabsolute==1), pdp=abs(out.ss.prel); else, pdp=out.ss.prel; end
    if(par.out.ss.plotabsolute==1), pdq=abs(out.ss.qrel); else, pdq=out.ss.qrel; end
    subplot(2,1,1), plot(out.ss.ttcom,pdp*100), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Pressure Relative Difference (%)','fontweight','bold'), else
        title('Pressure Relative Difference (%)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.ss.ttcom,pdq*100), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Flow Relative Difference (%)','fontweight','bold'), else
        title('Flow Relative Difference (%)','fontweight','bold'), end 
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\8rel.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\8rel.eps']), end
end
end

if(par.out.ss.plotcomps==1)
%compression ratios, discharge pressure setpoints, compressor power
    f9=figure(9); clf
    set(f9,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    tv=kron(out.ss.ccopt_int',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    legend(num2str([1:out.ss.n0.nc]'),'Location','SouthEast'), 
    title('Compression Ratios','FontWeight','bold'), 
    subaxis(1,3,2,'SpacingHoriz',0.05), 
    tv=kron(out.ss.csetopt_int',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    legend(num2str([1:out.ss.n0.nc]'),'Location','SouthEast'), 
    if(par.out.ss.units==1), title('Discharge Pressure Setpoints (psi)','FontWeight','bold'), else
        title('Discharge Pressure Setpoints (MPa)','FontWeight','bold'), end
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05), 
    tv=kron(out.ss.cpowopt_int'/1000,[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), 
    xlabel('hours'), ylabel(''), legend(num2str([1:out.ss.n0.nc]')), 
    title('Compressor Power (1000 hp)','FontWeight','bold'),
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\9comps.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\9comps.eps']), end
end

if(par.out.ss.plotdifferentials==1)
%compression ratios, discharge pressure setpoints, compressor power
    f10=figure(10); clf
    set(f10,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.ss.tt0,out.ss.ppoutopt-out.ss.ppinopt,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    legend(num2str([1:out.ss.n0.ne]'),'Location','SouthEast'), 
    if(par.out.ss.units==1), title('Pressure Differentials (psi)','FontWeight','bold'),  else
        title('Pressure Differentials (MPa)','FontWeight','bold'),  end
    subaxis(1,3,2,'SpacingHoriz',0.05), 
    plot(out.ss.tt0,out.ss.qqoutopt-out.ss.qqinopt,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    legend(num2str([1:out.ss.n0.ne]'),'Location','SouthEast'), 
    if(par.out.ss.units==1), title('Flow Differentials (mmscfd)','FontWeight','bold'), else
        title('Flow Differentials (kg/s)','FontWeight','bold'), end
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05),     
    tv=kron((out.ss.lmpout_int-out.ss.lmpin_int)',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), 
    xlabel('hours'), ylabel(''), legend(num2str([1:out.ss.n0.ne]')), 
    if(par.out.ss.units==1), title('Price Differentials ($/mscf)','FontWeight','bold'), else
        title('Price Differentials ($/mscf)','FontWeight','bold'), end
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\10differentials.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\10differentials.eps']), end
end

if(par.out.ss.plotcompshadow==1)
    %flexible gas withdrawals, lmps
    f11=figure(11); clf
    set(f11,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,2,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    tv=kron((out.ss.mult0_pmax_int.*(abs(out.ss.mult0_pmax_int)>10^(-3)))',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    legend(num2str([1:out.ss.n0.nc]'),'Location','SouthEast'), 
    if(par.out.ss.units==1), title('Shadow Price of Max Pressure ($/Psi/hr)','fontweight','bold'), else
    title('Shadow Price of Max Pressure ($/MPa/hr)','FontWeight','bold'), end
    subaxis(1,2,2,'MarginRight',0.05,'SpacingHoriz',0.05),
    tv=kron((out.ss.mult0_cmax_int.*(abs(out.ss.mult0_cmax_int)>10^(-3)))',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    title('Shadow Price of Max Comp. Power ($/hp)','fontweight','bold')
     legend(num2str([1:out.ss.n0.nc]'),'Location','SouthEast'), 
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\11compshadow.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\11compshadow.eps']), end
end

if(par.out.ss.plotaccuracy==1)
    f12=figure(12); clf
    set(f12,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,2,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.ss.flowbalrel'*100,'LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    legend(num2str([1:out.ss.n0.nv]'),'Location','SouthEast'), 
    title('Relative Nodal Flow Imbalance (%)','FontWeight','bold')
    subaxis(1,2,2,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.ss.flowbal','LineWidth',3), axis('tight'), xlabel('hours'), ylabel(''), 
    legend(num2str([1:out.ss.n0.nv]'),'Location','SouthEast'), 
    title('Absolute Nodal Flow Imbalance (mmscfd)','FontWeight','bold')
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\12accuracy.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\12accuracy.eps']), end
end

if(par.out.ss.plotpipemass==1)
    %optimization plot
    f13=figure(13);
    set(f13,'position',[20,400,800,400],'Color',[1 1 1]);
    tv=kron(out.ss.pipe_mass_int',[1 1])';
    plot(tp,tv,'LineWidth',3), axis('tight'), xlabel('time (hours)'), ylabel(''), 
    if(par.out.ss.units==1), title('Optimization Solution Mass in Pipe (mmscf)','fontweight','bold'), else
        title('Optimization Solution Mass in Pipe (kg)','fontweight','bold'), end   
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\5opt.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\13mass.eps']), end
end

if(par.out.ss.plotnetwork==1)
    f14=figure(14); clf
    plotpos=[20,100,800,800]; paperpos=[0 0 800 800]; papersize=[800 800];
    set(f14,'position',plotpos,'Color',[1 1 1]);
    gas_model_plotter_new(out.ss.n0); title('Original Network Topology')
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\13network.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\13network.eps']), end

    f15=figure(15); clf
    plotpos=[20,100,800,800]; paperpos=[0 0 800 800]; papersize=[800 800];
    set(f15,'position',plotpos,'Color',[1 1 1]);
    gas_model_plotter_new(out.ss.n);   title('Discretized Network Topology')
    if(par.out.ss.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\13network.pdf']), end
    if(par.out.ss.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\13network.eps']), end
end


%par.out=out;
%par.tr=tr;


%% backup
% %simulation pressure minus minimum
% f4=figure(4);
% set(f4,'position',[20,400,800,400],'Color',[1 1 1]);
% subplot(2,1,1), plot(tt/3600,pp-kron(ones(length(tt),1),par.n.p_min'/1000000)), 
% title('pressure minus minimum (MPa)')
% hold on, plot([tt(1)/3600 tt(end)/3600],[0 0])
% subplot(2,1,2), plot(tt/3600,qq), title('flow (kg/s)')
