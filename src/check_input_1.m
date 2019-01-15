function [par]=check_input_1(par)
%
% Anatoly Zlotnik, September 2017
% 
% Check parameter, options, and static network input
par.out.maxnv=100; par.out.maxne=100; par.out.maxnc=20; par.out.maxng=100;  %maximum nodes, edges, compressors, gnodes

if(par.out.intervals>0 && (par.out.intervals~=par.tr.optintervals ||...
        par.out.intervals~=par.out.intervals_out || par.tr.optintervals~=par.out.intervals_out))
    par.flag=1; par.message='Optimization intervals must equal I/O intervals if I/O intervals > 0';
return; end

if(par.ss.n0.nv>par.out.maxnv || par.ss.n0.ne>par.out.maxne ...
        || par.ss.n0.nc>par.out.maxnc || par.ss.n0.ng>par.out.maxng), 
    par.flag=1; par.message='System too large'; return; end

par.flag=0;
