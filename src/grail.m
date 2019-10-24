%% Gas Transient Optimization and Simulation
% Anatoly Zlotnik, updated September 2019
function [par]=grail(fpath)

% 1) steady state optimization
% 2) transient optimization
% 3) simulation
% 4) validation

% options:
% a) compressor power minimization
% b) economic optimization

% optimization solver: IPOPT
% ode solver: ode15i

%% Input data for problem definition
%clear

if(nargin==0)
    fnameid=fopen('model_folder.txt');
    fname = textscan(fnameid, '%s');
    fclose(fnameid);
    par.mfolder=fname{1}{1};
elseif(nargin>0)
    par.mfolder=fpath;
end


par=options_input(par);

%load static model (nodes, pipes, comps, gnodes) from xls fliles
if(exist([par.mfolder '\input_network.txt'])~=2 || par.out.update_from_xls==1)
    if(exist([par.mfolder '\input_network_nodes.xls'])==2), [par.ss.n0]=gas_model_reader_xls(par); else
    [par.ss.n0]=gas_model_reader_csv(par); end      
end

%first input check
par=check_input_1(par);
if(par.flag==1), disp(par.message), return; end



%% steady state solve
if(par.out.doss==1 || par.out.dosim==1)
                    
[par.ss.n0]=gas_model_reader_new(par.mfolder);          % load data from text file
if(par.ss.n0.nv>par.out.maxnv || par.ss.n0.ne>par.out.maxne || par.ss.n0.nc>par.out.maxng || par.ss.n0.ng>par.out.maxng), return; end

%if(par.ss.n0.nv>20), return; end;
[par.ss.n]=gas_model_reconstruct_new(par.ss.n0,par.ss.lmax,1);     % add extra nodes so segment lengths
                                                                   % are all below lmax
%f1=figure(1); clf
%gas_model_plotter_new(par.ss.n0);      % take a look at the network

%optimization parameters
par.ss.m.cdw=50;                      %comp. ratio derivative penalty scale
par.ss.m.ddw=50;                      %flex demand derivative penalty scale
par.ss.m.odw=100;                     %objective scale
par.ss.m.maxiter=400;                 %
par.ss.m.opt_tol=1e-6;                %

%model specifications
[par.ss]=model_spec(par.ss);

%demand function specifications
par.ss=econ_spec(par.ss,par.mfolder);

% Solve optimization
[par.ss]=tran_opt_base(par.ss);

%exit if not solved
if(par.ss.ip_info.status~=0), disp('Steady state optimization not feasible'), 
    fid=fopen([par.mfolder '\output_log.txt'],'w');
    fprintf(fid,['Steady-state solve status: ' num2str(par.ss.ip_info.status) '\n']);
    fclose(fid);
end

%process steady-state output
if(par.out.steadystateonly==1), [par]=process_output_ss(par); if(par.out.intervals_out>0), gas_out_plots_i(par); end; return; end
if(par.out.ss_check_exit==1), return; end
end
%% transient solve

[par.tr.n0]=gas_model_reader_new(par.mfolder);                   % load data from text file
if(par.ss.n0.nv>par.out.maxnv || par.ss.n0.ne>par.out.maxne || par.ss.n0.nc>par.out.maxng || par.ss.n0.ng>par.out.maxng), return; end

[par.tr.n]=gas_model_reconstruct_new(par.tr.n0,par.tr.lmax,0);     % add extra nodes so segment lengths
                                                                 % are all below lmax
%f1=figure(1); clf
%gas_model_plotter_new(par.tr.n0);      % take a look at the network

%optimization parameters
par.tr.m.cdw=50;                      %comp. ratio derivative penalty scale
par.tr.m.ddw=50;                      %flex demand derivative penalty scale
% par.tr.m.cdw=max(par.tr.Nvec);                      %comp. ratio derivative penalty scale
% par.tr.m.ddw=max(par.tr.Nvec);                      %flex demand derivative penalty scale

%model specifications
[par.tr]=model_spec(par.tr);

%demand function specifications
par.tr=econ_spec(par.tr,par.mfolder);


%% Solve optimization
[par.tr]=tran_opt_base(par.tr);

%% prepare simulation parameters
if(par.out.dosim==1)
    
%optimization solutions for setup
%load([par.mfolder '_' num2str(par.tr.lmax) 'km_' num2str(par.tr.Nvec(end)+1) 'p'],'par')
par.sim=par.ss;

par.sim.rtol0=1e-2; par.sim.atol0=1e-1;
par.sim.rtol1=1e-3; par.sim.atol1=1e-2;
par.sim.rtol=1e-5; par.sim.atol=1e-3;  %error tolerances for simulation
%par.sim.startup=1/4.2;     %startup time (fraction of horizon)
par.sim.startup=1/8;        %startup time (fraction of horizon)
par.sim.nperiods=3;         %number of periods after startup
par.sim.solsteps=96;        %solution steps per period
par.sim.fromss=1;

[par]=tran_sim_setup(par);

% execute simulation
[par.sim]=tran_sim_base(par.sim);

end
%% generate, save, and plot output

[par]=process_output_tr(par);
[par]=process_output_ss(par);
if(par.out.intervals_out==0), gas_out_plots(par); end
if(par.out.intervals_out>0), gas_out_plots_i(par); end

%pause(inf)
if(par.out.closeafter==1), close all, exit, end
