function [par]=tran_sim_setup(par)

sim=par.sim;

sim.m.Yd=sim.m.Yq;
sim.m.Yd(1:sim.m.fn,:)=par.tr.m.Yq(1:sim.m.fn,:);
sim.m.Ygd=interp1qr(par.tr.tt0,par.tr.fd0(1:length(par.tr.m.gd),:)',par.tr.m.xd')';
sim.m.Ygs=interp1qr(par.tr.tt0,par.tr.fd0(length(par.tr.m.gd)+1:length(par.tr.m.gall),:)',par.tr.m.xd')';
sim.m.Yf=interp1qr(par.tr.tt0,par.tr.fd0',par.tr.m.xd')';
sim.m.Yd(par.tr.m.guniqueind,:)=sim.m.Yd(par.tr.m.guniqueind,:)+...
    par.tr.m.gtod*sim.m.Yf;
sim.m.Yd1=mean(sim.m.Yd,2);
sim.m.xd=par.tr.m.xd;
sim.m.cc0=par.tr.cc0;
sim.m.Ys=par.tr.m.Pslack;
sim.m.N=par.tr.m.N;
sim.m.t=par.tr.m.t;
sim.m.x=par.tr.m.x;
sim.m.N1=par.tr.m.N1;
sim.m.tk=par.tr.m.tk;
sim.m.xk=par.tr.m.xk;

sim.tstep=sim.c.T/sim.solsteps;

sim.startupgrid=[0:sim.tstep/sim.c.Tsc:sim.m.Ts*sim.startup];
sim.periodgrid=[0:sim.tstep/sim.c.Tsc:sim.m.Ts];
sim.cyclesgrid=[0:sim.tstep/sim.c.Tsc:(sim.nperiods-1)*sim.m.Ts];
par.sim=sim;