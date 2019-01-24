function [par]=options_input(par)
%Anatoly Zlotnik, modified September 2017

%steady-state solve and simulation options
par.out.doss=1;     %do steady-state solve?
par.out.dosim=0;    %do simulations?
par.out.simstart=1; %0=start from transient IC, 1=start from ss solution

%read input parameter spreadsheet
if(exist([par.mfolder '\input_param.xls'])==2)
[inputpar]=xlsread([par.mfolder '\input_param.xls']);  %use if xls
else
inputpar_0=importdata([par.mfolder '\input_param.csv']);
inputpar(1:8,10)=inputpar_0.data(17:24);
inputpar(1:6,1)=inputpar_0.data(1:6);
inputpar(1:6,7)=inputpar_0.data(11:16);
inputpar(1:4,4)=inputpar_0.data(7:10);
end

par.tr.output_file=[par.mfolder '\print_log_ss.txt'];
par.ss.output_file=[par.mfolder '\print_log_tr.txt'];

%optimization options
if(inputpar(1)==0), par.tr.Nvec=0; end
if(inputpar(1)>1), 
par.tr.optintervals=inputpar(1,10);
par.tr.Nvec=[2.^([2:log2((par.tr.optintervals-1)/2)]) par.tr.optintervals-1]; end
par.tr.lmax=inputpar(2,10);     %artificial nodes are added so segments have length < lmax
par.tr.m.econweight=inputpar(3,10);      % value in [0,1] is weight on econ problem (rest on efficiency)
par.tr.m.maxiter=inputpar(4,10);         % max iterations
par.tr.m.opt_tol=10^inputpar(5,10);      % solver tolerance
par.tr.m.odw=10^inputpar(6,10);          % objective scale
par.tr.m.extension=inputpar(7,10);       % extend time horizon by this many optimization intervals
par.out.ss_check_exit=inputpar(8,10);    % exit after steady-state feasibility check

%physical system parameters
par.tr.c.gasT=inputpar(1,1);                      %gas temperature (K)
par.tr.c.gasG=inputpar(2,1);                     %gas gravity 
%gas gravity is relative to air = 28.9626. methane = 16.043, ethane = 30.070
universalR=8314.472;
par.tr.c.gasR=universalR/(28.9626*par.tr.c.gasG);
par.tr.c.a=sqrt(par.tr.c.gasR*par.tr.c.gasT);
par.tr.c.mpow=(inputpar(3,1)-1)/inputpar(3,1);  	%specific heat capacity ratio
par.tr.m.fuelfactor=inputpar(4,1);         %fuel factor for compressors
par.tr.c.T=3600*inputpar(5,1);                   %time horizon
par.out.doZ=inputpar(6,1);      %include compressibility?
par.tr.m.doZ=par.out.doZ;

%input options
par.out.units=inputpar(1,4);     %for I/0 and display. 0 = SI, 1 = standard
%plot pressure in psi if units = 1, MPa if 0
%plot pressure in mmscfd if units = 1, kg/s if 0
par.tr.units=par.out.units; par.ss.units=par.out.units;
par.out.intervals=inputpar(2,4);  %input intervals
par.tr.intervals=par.out.intervals; par.ss.intervals=par.out.intervals;
par.tr.m.intervals=par.tr.intervals;
par.out.update_from_xls=inputpar(3,4);     %update stadic network. 0 = no, 1 = yes
par.out.use_init_state=inputpar(4,4);     %use provided initial state. 0 = no, 1= yes
par.tr.m.use_init_state=par.out.use_init_state;


%output options
par.out.savecsvoutput=inputpar(1,7);
par.out.steadystateonly=inputpar(2,7);
par.out.plotmarketflowlims=inputpar(3,7);
par.out.plotmarketpricebids=inputpar(3,7);
par.out.plotsolflow=inputpar(3,7);
par.out.plotsolprice=inputpar(3,7);
par.out.plotopt=inputpar(3,7);
par.out.plotcomps=inputpar(3,7);
par.out.plotdifferentials=inputpar(3,7);
par.out.plotcompshadow=inputpar(3,7);
par.out.plotnetwork=inputpar(3,7);
par.out.plotpipemass=inputpar(3,7);
par.out.plotpdf=inputpar(4,7);
par.out.closeafter=inputpar(5,7);
par.out.plotaccuracy=0;
par.out.intervals_out=inputpar(6,7);     %output intervals
%par.out.plotaccuracy=inputpar(22);
%par.out.plotnetwork=inputpar(23);
%par.out.plotpdf=inputpar(24);
%par.out.closeafter=inputpar(25);

par.out.plotlarge=0;
if(par.out.dosim==1)
par.out.plotsim=1;
par.out.plotcompabs=1;
par.out.plotcomprel=1;
end
par.out.ploteps=0;


%plotting options
par.out.plotnodal=1;            %only plot nodal (n0) if 1, plot all (n) if 9
par.out.plotabsolute=0;         %plot absolute distance |x - y|

%problem type
par.tr.m.pdcstr=1;              % exact discharge pressure constraint
par.tr.m.cpowcstr=1;            % compressor power constraint

%state save points
if(exist([par.mfolder '\input_state_save_pts.csv'])==2)
    par.tr.m.save_state=1; par.ss.m.save_state=1;
    par.tr.m.state_save_pts=double(dlmread([par.mfolder '\input_state_save_pts.csv']));
elseif(exist([par.mfolder '\input_state_save_pts.csv'])==0)
    par.tr.m.save_state=0; par.ss.m.save_state=0;
end

%time discretization
if(par.out.doss==1 || par.out.dosim==1)
    par.ss.Nvec=[0];    % vector of time discretization orders for optimization sequence (steady state)
    par.ss.lmax=1;      % max segment length (km) (steady state)
    if(par.out.dosim==0) par.ss.lmax=1000; end
    par.ss.c.T=par.tr.c.T;
    par.ss.c.gasR=par.tr.c.gasR;                      %gas constant R (J/kg K)
    par.ss.c.gasT=par.tr.c.gasT;                      %gas temperature (K)
    par.ss.c.gasG=par.tr.c.gasG;                     %gas gravity 
    par.ss.c.a=par.tr.c.a; par.ss.c.mpow=par.tr.c.mpow; par.ss.m.fuelfactor=par.tr.m.fuelfactor;     
    par.ss.m.econweight=par.tr.m.econweight;
    par.ss.m.pdcstr=par.tr.m.pdcstr; par.ss.m.cpowcstr=par.tr.m.cpowcstr;
    par.ss.m.doZ=par.out.doZ;
    par.ss.m.use_init_state=0;
    par.ss.m.state_save_pts=1;
    par.ss.m.extension=0;
end
%par.tr.m.compints=6;  % intervals in the period with constant compression

%physical parameters
