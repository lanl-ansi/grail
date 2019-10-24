function [par]=model_spec(par)
%Anatoly Zlotnik, modified February 2017

%rescaling parameters
par.c.psi_to_pascal=6894.75729;
par.c.km_to_m=1000;
universalR=8314.472;    %J/k-mol * K
molecwghtAir=28.9626;   %Molecular weight of air (g mol)
oneATM=101.325;     %Pascal
%par.c.mmscfd_to_kgps=0.24085926;   convert at given gas gravity and temp. below
par.c.mmscfd_to_kgps=(10^3)*0.02832/86400*(oneATM/(universalR*par.c.gasT)*(100)^3)*par.c.gasG*molecwghtAir; %conversion notes
par.c.mmscfd_to_hp=1000/0.0593807;  %fix this later (should depend also on calorific value)
%par.c.inch_to_m=2.54/100;
%par.c.miles_to_m=1609.34;


%CNGA parameters ( Z=1/(b1+b2*p) ) in with metric I/O
a1=344400; a2=1.785; a3=3.825;
b1=1+a1*(101.325*1000/par.c.psi_to_pascal)*10^(a2*par.c.gasG)/(1.8*par.c.gasT)^a3;
b2=a1*10^(a2*par.c.gasG)/par.c.psi_to_pascal/(1.8*par.c.gasT)^a3;
if(par.m.doZ==1), par.c.b1=b1; par.c.b2=b2;
else par.c.a=sqrt(par.c.gasR*par.c.gasT); end
%else par.c.a=sqrt(par.c.gasR*par.c.gasT/(b1+b2*mean([par.n.p_min;par.n.p_max]))); end


%if(~exist('par.m'))
%change units: km -> m
par.n.pipe_length=par.n.pipe_length*par.c.km_to_m;
%end

%additional mpow definition
par.m.mpow=par.c.mpow;             

%nondimensionalization
par.c.psc=par.n.p_min(par.n.slack_nodes(1));    %nominal pressure
par.c.R=1000;                                   %nondim. length (m)
par.c.Tsc=par.c.R/par.c.a;                      %time rescaling
par.c.Xsc=par.c.R;                              %space rescaling
par.c.qsc=par.c.psc/par.c.a;                    %flow rescaling

%rescaled parameters
%if(~exist('par.c.T')), par.c.T=86400; end
par.m.Ts=par.c.T/par.c.Tsc;               %rescaled nondimensionalized time horizon
par.m.p_min_nd=par.n.p_min/par.c.psc;     %nondimensionalized min pressure
par.m.p_max_nd=par.n.p_max/par.c.psc;     %nondimensionalized min pressure

%construct model parameter data
par.m.Lk=sparse(par.n.pipe_length);         %pipe lengths
par.m.Dk=sparse(par.n.diameter);            %pipe diameter
par.m.lamk=sparse(par.n.lambda);            %pipe friction factor
par.m.ML=diag(sparse(par.m.Lk/4/par.c.R));  %parameter matrix "Lambda"
par.m.Ma=diag(sparse(-par.c.R./par.m.Lk));  %parameter matrix "K"
par.m.Xs=diag(sparse(par.m.Dk.^2/4*pi));    %cross section area matrix
par.m.xs=sparse(par.m.Dk.^2/4*pi);      %cross section area vector
par.m.FN=length(par.n.nonslack_nodes);  % # of flow nodes
par.m.PN=length(par.n.slack_nodes);     % # of pressure nodes
par.m.NE=par.n.ne;                      % # of edges
par.m.M=par.m.FN+par.m.NE;              % # of state variables - edges and nodes
par.m.C=size(par.n.loc_node,1);         % # of compressors
par.m.dnodes=par.n.nonslack_nodes;      % demand nodes
par.m.snodes=par.n.slack_nodes;         % demand nodes
par.m.comp_pos=par.n.comp_pos;          % compressor loc_node and to_edge
par.m.Am=par.n.Am;                      % full adjacency matrix
par.m.Ad=par.n.Ad;                      % demand adjacency submatrix
par.m.Adp=(par.m.Ad>0).*par.m.Ad;       % positive demand adjacency submatrix
par.m.Adn=(par.m.Ad<0).*par.m.Ad;       % negative demand adjacency submatrix
par.m.Amp=(par.m.Am>0).*par.m.Am;       % positive demand adjacency submatrix
par.m.Amn=(par.m.Am<0).*par.m.Am;       % negative demand adjacency submatrix
par.m.As=par.n.As;                      % supply adjacency submatrix
par.m.R=par.c.R;                        % nondimensionalization length
if(par.m.doZ==1)
par.m.b1=par.c.b1; par.m.b2=par.c.b2;   % CNGA parameters
par.m.psc=par.c.psc;                    % pressure rescaling
end

%determine if compressors are at demand or slack nodes
par.m.spos=[];              %compressors at slack nodes
par.m.dpos=[];              %compressors at demand nodes
for j=1:par.m.C
    if(par.n.isslack(par.n.loc_node(j))==1)
        par.m.spos=[par.m.spos;j];                  %slackbus compressor positions
    else
        par.m.dpos=[par.m.dpos;j];                  %non slackbus compressor positions
    end
end

%determine compressor location index in slack nodes enumeration
dvec=zeros(par.n.nv,1); dvec(par.m.comp_pos(:,1))=[1:par.m.C]'; dvec(par.m.dnodes)=[];
par.m.sppos=[];
for j=1:length(par.m.spos)
    par.m.sppos(j,:)=find(dvec==par.m.spos(j));
end

%determine compressor location index in demand nodes enumeration
dvec=zeros(par.n.nv,1); dvec(par.m.comp_pos(:,1))=[1:par.m.C]'; dvec(par.m.snodes)=[];
par.m.ppos=[];
for j=1:length(par.m.dpos)
    par.m.ppos(j,:)=find(dvec==par.m.dpos(j));
end


% %compressor feasibility constraints
% par.m.eff=par.m.fuelfactor/(max(par.n.c_max)^par.m.mpow-1); 	%fuel factor
% %par.m.eff=0.00075; 	%fuel factor
% %par.m.eff=0.00025;
% if(par.m.doZ==1), Zsc=(par.c.b1+par.c.b2*par.n.p_min(par.n.comp_pos(:,1)))./(par.c.b1+par.c.b2*par.n.p_max(par.n.comp_pos(:,1))); else
%     Zsc=1; end
% par.m.boost_pow_max_nd=(1/par.m.eff)*par.n.hp_max/par.c.mmscfd_to_hp*par.c.mmscfd_to_kgps/par.c.qsc...
%     .*Zsc;

%compressor feasibility constraints
if(par.m.doZ==1), Zsc=(par.c.b1+par.c.b2*par.n.p_min(par.n.comp_pos(:,1)))./(par.c.b1+par.c.b2*par.n.p_max(par.n.comp_pos(:,1))); else
    Zsc=1; end
par.m.eff=par.m.fuelfactor./((par.n.c_max.*Zsc).^par.m.mpow-1); 	%fuel factor

%par.m.boost_pow_max_nd=(1./par.m.eff).*par.n.hp_max/par.c.mmscfd_to_hp*par.c.mmscfd_to_kgps/par.c.qsc;
par.m.Wc=286.76/par.c.gasG*par.c.gasT*(par.c.gamm/(par.c.gamm-1));
par.m.boost_pow_max_nd=par.n.hp_max/par.m.Wc/par.c.qsc;

par.m.flow_min_nd=par.n.flow_min/par.c.qsc;
par.m.flow_max_nd=par.n.flow_max/par.c.qsc;
