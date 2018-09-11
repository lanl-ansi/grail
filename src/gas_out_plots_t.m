function []=gas_out_plots_t(par)

out=par.out; mfolder=out.mfolder;

if(par.out.plotlarge==1), 
    plotpos=[20,100,1600,860]; paperpos=[0 0 1600 860]; papersize=[1600 860]; else
    plotpos=[20,400,1200,400]; paperpos=[0 0 1200 400]; papersize=[1200 400]; end

int_bounds=out.int_bounds;

if(par.out.plotmarketflowlims==1)
    %fixed gas withdrawals, bounds on flexible demands and supplies
    f1=figure(1); clf
    set(f1,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    ts=timeseries([out.dbase_int(:,out.guniqueind);out.dbase_int(end,out.guniqueind)]); ts=setinterpmethod(ts, 'zoh');
    ts.Time=int_bounds';
    plot(ts,'LineWidth',3), axis('tight'), xlabel('hours')
    if(par.out.units==1), title('Baseline Gas Withdrawals (mmscfd)','fontweight','bold'), else
        title('Baseline Gas Withdrawals (kg/s)','fontweight','bold'), end
    if(length(out.guniqueind)<20), legend(num2str(out.gunique),'Location','SouthEast'), end
    subaxis(1,3,2,'SpacingHoriz',0.05), 
    plot(out.td,out.gdub,'LineWidth',3), axis('tight'), xlabel('hours'), hold on,
    plot(out.td,out.gdlb,'--','LineWidth',3), hold off, legend(num2str(out.gd),'Location','SouthEast')
    if(par.out.units==1), title('Bounds on Buyer Offtakes (mmscfd)','fontweight','bold'), else
    title('Bounds on Buyer Offtakes (kg/s))','FontWeight','bold'), end
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05), 
    plot(out.td,out.gsub,'LineWidth',3), axis('tight'), xlabel('hours'), hold on,
    plot(out.td,out.gslb,'--','LineWidth',3), hold off, legend(num2str(out.gd),'Location','SouthEast')
    if(par.out.units==1), title('Bounds on Seller Supplies (mmscfd)','fontweight','bold'), else
    title('Bounds on Seller Supplies (kg/s))','FontWeight','bold'), end
    if(par.out.plotpdf==1)
        set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
            'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
        eval(['print -dpdf ' mfolder '\1marketflowlims.pdf']), end
    if(par.out.ploteps==1)
        set(gcf,'PaperPositionMode','auto')
        eval(['print -depsc ' mfolder '\1marketflowlims.eps'])
    end
end

if(par.out.plotmarketpricebids==1)
    %price bids and offers
    f2=figure(2); clf
    set(f2,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
	plot(out.tt0,out.Prslack,'LineWidth',3), axis('tight'), xlabel('hours')
    if(par.out.units==1), title('Pressure (slack) Node Prices ($/mscf)','fontweight','bold'), else
    title('Pressure (slack) Node (10$/kg)','FontWeight','bold'), end
    legend(num2str(out.PN),'Location','SouthEast')
    subaxis(1,3,2,'SpacingHoriz',0.05), plot(out.td,out.Prd,'LineWidth',3), axis('tight'), xlabel('hours')
    if(par.out.units==1), title('Demand Bid Prices ($/mscf)','fontweight','bold'), else
    title('Demand Bid Prices (10$/kg)','FontWeight','bold'), end
    legend(num2str(out.gd),'Location','SouthEast')
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05), 
    plot(out.td,out.Prs,'LineWidth',3), axis('tight'), xlabel('hours')
    if(par.out.units==1), title('Supply Offer Prices ($/mscf)','fontweight','bold'), else
    title('Supply Offer Prices (10$/kg)','FontWeight','bold'), end
    legend(num2str(out.gs),'Location','SouthEast')
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\2marketpricebids.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\2marketpricebids.eps']), end
end


if(par.out.plotsolflow==1)
    %flow schedule solution
    f3=figure(3); clf
    set(f3,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.tt0,[-out.supp_flow out.dgflows],'LineWidth',3), xlabel('hours'), hold off
    if(par.out.units==1), title('Cleared Nodal Gas Withdrawals (mmscfd)','fontweight','bold'), else
        title('Cleared Nodal Gas Withdrawals (kg/s)','fontweight','bold'), end
    legend(num2str([out.PN;out.gunique]),'Location','SouthEast')
    subaxis(1,3,2,'SpacingHoriz',0.05),
    plot(out.tt0,out.gdsol,'LineWidth',3), axis('tight'), xlabel('hours'),
    if(par.out.units==1), title('Demand gNode Purchases (mmscfd)','fontweight','bold'), else
        title('Demand gNode Purchases  (kg/s)','fontweight','bold'), end
    legend(num2str(out.gd),'Location','SouthEast')
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05),
    plot(out.tt0,[out.supp_flow out.gssol],'LineWidth',3), axis('tight'), xlabel('hours'),
    if(par.out.units==1), title('Supply gNode Sales (mmscfd)','fontweight','bold'), else
        title('Supply gNode Sales  (kg/s)','fontweight','bold'), end
    legend(num2str([out.PN;out.gs]),'Location','SouthEast')
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\3solflow.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\3solflow.eps']), end
end

if(par.out.plotsolprice==1)
    %flexible gas withdrawals, lmps
    f4=figure(4); clf
    set(f4,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,2,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.tt0,out.trlmp,'LineWidth',3), axis('tight'), xlabel('hours')
    if(par.out.units==1), title('All LMPs ($/mscf)','fontweight','bold'), else
    title('All LMPs (10$/kg)','FontWeight','bold'), end
    subaxis(1,2,2,'MarginRight',0.05,'SpacingHoriz',0.05),
    plot(out.tt0,out.dglmp,'LineWidth',3), axis('tight'), xlabel('hours')
    if(par.out.units==1), title('gNode LMPs ($/mscf)','fontweight','bold'), else
    title('gNode LMPs (10$/kg)','FontWeight','bold'), end
    legend(num2str(out.gunique),'Location','SouthEast')
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\4solprice.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\4solprice.eps']), end
end

if(par.out.plotopt==1)
    %optimization plot
    f5=figure(5);
    set(f5,'position',[20,400,800,400],'Color',[1 1 1]);
    subplot(2,1,1), plot(out.tt0,out.ppopt), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Optimization Solution Nodal Pressures (psi)','fontweight','bold'), else
        title('Optimization Solution Nodal Pressures (MPa)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.tt0,out.qqopt), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Optimization Solution Flows (mmscfd)','fontweight','bold'), else
        title('Optimization Solution Flows (kg/s)','fontweight','bold'), end   
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\5opt.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\5opt.eps']), end
end

if(par.out.dosim==1)
if(par.out.plotsim==1)
    %simulation plot
    f6=figure(6);
    set(f6,'position',[20,400,800,400],'Color',[1 1 1]);
    subplot(2,1,1), plot(out.tt,out.ppsim), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Simulation Solution Nodal Pressures (psi)','fontweight','bold'), else
        title('Simulation Solution Nodal Pressures (MPa)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.tt,out.qqsim), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Simulation Solution Flows (mmscfd)','fontweight','bold'), else
        title('Simulation Solution Flows (kg/s)','fontweight','bold'), end
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\6sim.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\6sim.eps']), end
end

if(par.out.plotcompabs==1)
    %comparison plot absolute
    f7=figure(7);
    set(f7,'position',[20,400,800,400],'Color',[1 1 1]);
    if(par.out.plotabsolute==1), pdp=abs(out.pdiff); else, pdp=out.pdiff; end
    if(par.out.plotabsolute==1), pdq=abs(out.qdiff); else, pdq=out.qdiff; end
    subplot(2,1,1), plot(out.ttcom,pdp), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Pressure Absolute Difference (psi)','fontweight','bold'), else
        title('Pressure Absolute Difference (MPa)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.ttcom,pdq), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Flow Absolute Difference (mmscfd)','fontweight','bold'), else
        title('Flow Absolute Difference (kg/s)','fontweight','bold'), end   
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\7diff.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\7diff.eps']), end
end

if(par.out.plotcomprel==1)
    %comparison plot relative
    f8=figure(8);
    set(f8,'position',[20,400,800,400],'Color',[1 1 1]);
    if(par.out.plotabsolute==1), pdp=abs(out.prel); else, pdp=out.prel; end
    if(par.out.plotabsolute==1), pdq=abs(out.qrel); else, pdq=out.qrel; end
    subplot(2,1,1), plot(out.ttcom,pdp*100), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Pressure Relative Difference (%)','fontweight','bold'), else
        title('Pressure Relative Difference (%)','fontweight','bold'), end   
    subplot(2,1,2), plot(out.ttcom,pdq*100), axis('tight'), xlabel('time (hours)')
    if(par.out.units==1), title('Flow Relative Difference (%)','fontweight','bold'), else
        title('Flow Relative Difference (%)','fontweight','bold'), end 
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\8rel.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\8rel.eps']), end
end
end

if(par.out.plotcomps==1)
%compression ratios, discharge pressure setpoints, compressor power
    f9=figure(9); clf
    set(f9,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.tt0,out.ccopt,'LineWidth',3), axis('tight'), xlabel('hours')
    legend(num2str([1:out.n0.nc]'),'Location','SouthEast'), 
    title('Compression Ratios','FontWeight','bold'), 
    subaxis(1,3,2,'SpacingHoriz',0.05), 
    plot(out.tt0,out.csetopt,'LineWidth',3), axis('tight'), xlabel('hours')
    legend(num2str([1:out.n0.nc]'),'Location','SouthEast'), 
    if(par.out.units==1), title('Discharge Pressure Setpoints (psi)','FontWeight','bold'), else
        title('Discharge Pressure Setpoints (MPa)','FontWeight','bold'), end
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05), 
    plot(out.tt0,out.cpowopt/1000,'LineWidth',3), axis('tight'), 
    xlabel('hours'), legend(num2str([1:out.n0.nc]')), 
    title('Compressor Power (1000 hp)','FontWeight','bold'),
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\9comps.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\9comps.eps']), end
end

if(par.out.plotdifferentials==1)
%compression ratios, discharge pressure setpoints, compressor power
    f10=figure(10); clf
    set(f10,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,3,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.tt0,out.ppoutopt-out.ppinopt,'LineWidth',3), axis('tight'), xlabel('hours')
    legend(num2str([1:out.n0.ne]'),'Location','SouthEast'), 
    if(par.out.units==1), title('Pressure Differentials (psi)','FontWeight','bold'),  else
        title('Pressure Differentials (MPa)','FontWeight','bold'),  end
    subaxis(1,3,2,'SpacingHoriz',0.05), 
    plot(out.tt0,out.qqoutopt-out.qqinopt,'LineWidth',3), axis('tight'), xlabel('hours')
    legend(num2str([1:out.n0.ne]'),'Location','SouthEast'), 
    if(par.out.units==1), title('Flow Differentials (mmscfd)','FontWeight','bold'), else
        title('Flow Differentials (kg/s)','FontWeight','bold'), end
    subaxis(1,3,3,'MarginRight',0.05,'SpacingHoriz',0.05), 
    plot(out.tt0,out.lmpout-out.lmpin,'LineWidth',3), axis('tight'), 
    xlabel('hours'), legend(num2str([1:out.n0.ne]')), 
    if(par.out.units==1), title('Price Differentials ($/mscf)','FontWeight','bold'), else
        title('Price Differentials ($/mscf)','FontWeight','bold'), end
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\10differentials.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\10differentials.eps']), end
end

if(par.out.plotcompshadow==1)
    %flexible gas withdrawals, lmps
    f11=figure(11); clf
    set(f11,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,2,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.tt0,out.mult0_pmax,'LineWidth',3), axis('tight'), xlabel('hours')
    legend(num2str([1:out.n0.nc]'),'Location','SouthEast'), 
    if(par.out.units==1), title('Shadow Price of Max Pressure ($/Psi/hr)','fontweight','bold'), else
    title('Shadow Price of Max Pressure ($/MPa/hr)','FontWeight','bold'), end
    subaxis(1,2,2,'MarginRight',0.05,'SpacingHoriz',0.05),
    plot(out.tt0,out.mult0_cmax,'LineWidth',3), axis('tight'), xlabel('hours')
    title('Shadow Price of Max Comp. Power ($/hp)','fontweight','bold')
     legend(num2str([1:out.n0.nc]'),'Location','SouthEast'), 
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\11compshadow.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\11compshadow.eps']), end
end

if(par.out.plotaccuracy==1)
    f12=figure(12); clf
    set(f12,'position',plotpos,'Color',[1 1 1]);
    subaxis(1,2,1,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.flowbalrel'*100,'LineWidth',3), axis('tight'), xlabel('hours')
    legend(num2str([1:out.n0.nv]'),'Location','SouthEast'), 
    title('Relative Nodal Flow Imbalance (%)','FontWeight','bold')
    subaxis(1,2,2,'MarginLeft',0.05,'SpacingHoriz',0.05), 
    plot(out.flowbal','LineWidth',3), axis('tight'), xlabel('hours')
    legend(num2str([1:out.n0.nv]'),'Location','SouthEast'), 
    title('Absolute Nodal Flow Imbalance (mmscfd)','FontWeight','bold')
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\12accuracy.pdf']), end
    if(par.out.ploteps==1)
    set(gcf,'PaperPositionMode','auto')
    eval(['print -depsc ' mfolder '\12accuracy.eps']), end
end

if(par.out.plotnetwork==1)
    f13=figure(13); clf
    plotpos=[20,100,800,800]; paperpos=[0 0 800 800]; papersize=[800 800];
    set(f13,'position',plotpos,'Color',[1 1 1]);
    gas_model_plotter_new(out.n0);
    if(par.out.plotpdf==1)
    set(gcf,'PaperPositionMode', 'manual','PaperUnits','points', ...
        'Paperposition',paperpos), set(gcf, 'PaperSize', papersize)
    eval(['print -dpdf ' mfolder '\13network.pdf']), end
    if(par.out.ploteps==1)
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
