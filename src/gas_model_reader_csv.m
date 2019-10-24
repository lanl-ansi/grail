function [network]=gas_model_reader_csv(par)
% Anatoly Zlotnik, modified September 2019
%
% usage: [network]=gas_model_reader_xls(mod_num)
% function: read static network data from xls files and translate to txt format
%
%
% input: mfolder = name of folder containing network.txt
%
% input file format: csv files
% input_network_nodes: [node id] [node name] [x coord] [y coord] [pmin] [pmax] [qmin] [qmax] [slack bool]
% input_network_pipes: [pipe id] [pipe name] [from id] [to id] [diameter] [length] [friction factor] [disc seg]
% input_network_comps: [comp id] [comp name] [from id] [to id] [cmin] [cmax] [hp max] [min flow] [max flow]
% input_network_gnodes: [gnode id] [gnode name] [node id]
% note: numbering begins from 1
%
% output: network struct
%         network.nv = number of vertices
%         network.ne = number of edges
%         network.nc = number of compressors
% node data:
%         network.node_name = node name strings
%         network.xcoord = x coordinate
%         network.ycoord = y coordinate
%         network.p_min = min pressure at vertices (MPa \ psi)
%         network.p_max = max pressure at vertices (MPa \ psi)
%         network.q_min = min injection (negative = consumption) at vertices (kg/sec \ mmscfd (1 atm + 60F))
%         network.q_max = max injection (negative = consumption) at vertices (kg/sec \ mmscfd (1 atm + 60F))
%         network.isslack = 1 -> slack node, 0 -> not slack
%         network.nomflow = nominal flow, if 9th column in file
% edge data:
%         network.pipe_name = pipe name strings
%         network.from_id = start node #
%         network.to_id = end node #
%         network.diameter = pipe diameter (meters \ inches)
%         network.pipe_length = pipe length (km \ miles)
%         network.lambda = multiplier
%         network.disc_seg = segments for discretization by pipe (0 = lmax)
% compressor data:
%         network.comp_name = comp name strings
%         network.loc_node = compressor location node #
%         network.to_edge = edge # to which flow is boosted 
%         network.c_min = min compression ratio
%         network.c_max = max compression ratio     
%         network.hp_max = max compressor power (hp)
%         network.flow_min = min flow through compressor (kg/s \ mmscfd)
%         network.flow_max = max flow through compressor (kg/s \ mmscfd)
% gnode data:
%         network.gnode_name = gnode name strings
%         network.gnode_loc = gnode physical location

mfolder=par.mfolder;

%conversion factors
inch_to_m=0.0254;
miles_to_km=1.60934;
psi_to_pascal=6894.75729;
m_to_ft=3.28084;
hp_to_watt=745.7;
%mmscfd to kgps - see conversion notes.  
%this depends on gas composition and definition of "standard"
universalR=8314.472;    %J/k-mol * K
standardP=101.325;      %Pascal (=14.69psi)
standardT=288.706;      %Kelvin (=60F)
molecwghtAir=28.9626;   %Molecular weight of air (g mol)
mmscfd_to_kgps=(10^3)*(1/m_to_ft)^3/86400*(standardP/(universalR*standardT)*(100)^3)*par.tr.c.gasG*molecwghtAir;

%import from xls
%[node_data,node_text]=xlsread([mfolder '\input_network_nodes.xls']);
%[pipe_data,pipe_text]=xlsread([mfolder '\input_network_pipes.xls']);
%[comp_data,comp_text]=xlsread([mfolder '\input_network_comps.xls']);
%[gnode_data,gnode_text]=xlsread([mfolder '\input_network_gnodes.xls']);

%import from csv
nodedat=importdata([par.mfolder '\input_network_nodes.csv']);
pipedat=importdata([par.mfolder '\input_network_pipes.csv']);
compdat=importdata([par.mfolder '\input_network_comps.csv']);
gnodedat=importdata([par.mfolder '\input_network_gnodes.csv']);
node_text=nodedat.textdata; pipe_text=pipedat.textdata;
comp_text=compdat.textdata; gnode_text=gnodedat.textdata;
node_data=[zeros(size(nodedat.data,1),2) nodedat.data];
pipe_data=[zeros(size(pipedat.data,1),2) pipedat.data];
comp_data=[zeros(size(compdat.data,1),2) compdat.data];
gnode_data=[zeros(size(gnodedat.data,1),2) gnodedat.data];

nv=size(node_data,1);                   %number of vertices
ne=size(pipe_data,1)+size(comp_data,1); %number of edges
nc=size(comp_data,1);                   %number of compressors
ng=size(gnode_data,1);                  %number of gnodes

node_name=node_text(2:end,2);
pipe_name=[pipe_text(2:end,2);comp_text(2:end,2)];
comp_name=comp_text(2:end,2);
gnode_name=gnode_text(2:end,2);

%compressor pipe standard
if(par.out.units==0), stdL=.25; stdD=1;           %250 m and 1 m
else, stdL=.25/miles_to_km; stdD=1/inch_to_m;   %.15 miles and 39.4"
end
stdLambda=0.001;

%information for each node
xcoord=node_data(:,3);    %longitude
ycoord=node_data(:,4);    %lattitude
p_min=node_data(:,5);     %min pressure
p_max=node_data(:,6);     %max pressure
q_min=node_data(:,7);     %min withdrawal
q_max=node_data(:,8);     %max withdrawal
isslack=node_data(:,9);   %1 if slack bus otherwise zero

%information for each edge
from_id=[pipe_data(:,3); comp_data(:,3)];      %from node
to_id=[pipe_data(:,4); comp_data(:,4)];        %to node
diameter=[pipe_data(:,5); stdD*ones(nc,1)];    %diameter (m)
pipe_length=[pipe_data(:,6); stdL*ones(nc,1)]; %length (km)
lambda=[pipe_data(:,7); stdLambda*ones(nc,1)]; %friction coefficient 
disc_seg=[pipe_data(:,8); ones(nc,1)];         %additional discretization points by pipe

%information for each compressor
loc_node=comp_data(:,3);    %location node for compressor
to_edge=[size(pipe_data,1)+1:size(pipe_data,1)+nc]';       %to edge
c_min=comp_data(:,5);       %min compression ratio
c_max=comp_data(:,6);       %max compression ratio
hp_max=comp_data(:,7);         %max compressor power (MW)
flow_min=comp_data(:,8);         %min compressor flow (kg/s)
flow_max=comp_data(:,9);         %max compressor flow (kg/s)
comp_bool=[zeros(ne,1);ones(nc,1)];      % = 1 if compressor on the edge, 0 otherwise
comp_pos=[loc_node to_edge];

%information for each gnode
phys_node=gnode_data(:,3);  %physical node location of gnode

%adjacency matrix
Amm=sparse(nv,ne); Amp=sparse(nv,ne);
for i=1:ne
    Amm(from_id(i),i)=-1; Amp(to_id(i),i)=1;
end
Am=Amp+Amm;
% % for i=1:nv
% %    Ad(i,:)=-(from_id==i)+(to_id==i);
% % end
%slack and demand submatrices
%As=Am(isslack==1,:);
%Ad=Am(isslack==0,:);

slack_nodes=find(isslack);      %slack nodes
nonslack_nodes=setdiff([1:nv]',slack_nodes); %nodes that are not the slack bus

network.nv=nv;              %number of vertices
network.ne=ne;              %number of edges
network.nc=nc;              %number of compressors
network.ng=ng;              %number of compressors
network.node_name=node_name;      %node name strings
network.pipe_name=pipe_name;      %node name strings
network.comp_name=comp_name;      %node name strings
network.gnode_name=gnode_name;    %node name strings
network.xcoord=xcoord;      %x coordinate of vertices
network.ycoord=ycoord;      %y coordinate of vertices
network.p_min=p_min;        %min pressure at vertices (psi)
network.p_max=p_max;        %max pressure at vertices (psi)
network.isslack=isslack;    % = 1 if slack node
network.from_id=from_id;    %start node #
network.to_id=to_id;        %end node #
network.diameter=diameter;  %pipe diameter (m)
network.pipe_length=pipe_length;      %pipe length (km)
network.lambda=lambda;      %friction factor
network.disc_seg=disc_seg;  %additional discretization points by pipe
network.c_min=c_min;        %min compression ratio
network.c_max=c_max;        %max compression ratio
network.loc_node=loc_node;  %compressor location node #
network.to_edge=to_edge;    %edge # to which flow is boosted 
network.c_min=c_min;        %min compression ratio
network.c_max=c_max;        %max compression ratio     
network.hp_max=hp_max;      %max compressor power (hp)
network.flow_min=flow_min;  %min compressor flow (mmscfd)
network.flow_max=flow_max;  %max compressor flow (mmscfd)
network.comp_pos=comp_pos;  %compressor positions (to sub into adjacency matrix)
network.comp_bool=comp_bool;% = 1 if compressor on edge, = 0 otherwise
network.Am=Am;              %adjacency matrix
network.Amp=Amp;            %adjacency matrix positive part
network.Amm=Amm;            %adjacency matrix negative part
network.slack_nodes=slack_nodes; %slack nodes
network.nonslack_nodes=nonslack_nodes; %nonslack nodes
network.phys_node=phys_node;    %physical node location of gnode

% output file format: space separated file
% first row: [# nodes] [# edges] [# compressors] [#gnodes]
% node rows: [node #] [node name] [xcoord] [ycoord] [pmin] [pmax] [qmin] [qmax] [slack bool]
% edge rows: [edge #] [pipe name] [from node] [to node] [diameter (m)] [length (km)] [friction factor]
% comp rows: [comp #] [comp name] [loc node] [to edge] [cmin] [cmax] [hp max] [min flow] [max flow]
% gnode rows: [gnode #] [gnode name] [phys loc]

%converted units
if(par.out.units==1)    %if given in standard, convert to SI at given gas gravity
    p_min=p_min*psi_to_pascal; p_max=p_max*psi_to_pascal;
    q_min=q_min*mmscfd_to_kgps; q_max=q_max*mmscfd_to_kgps; 
    diameter=diameter*inch_to_m; pipe_length=pipe_length*miles_to_km;
    flow_min=flow_min*mmscfd_to_kgps;
    flow_max=flow_max*mmscfd_to_kgps;
    hp_max=hp_max*hp_to_watt;
end

fid = fopen([mfolder '\input_network.txt'],'w');

line1=[num2str(nv) ' ' num2str(ne) ' ' num2str(nc) ' ' num2str(ng) '\n'];
fprintf(fid,line1);

for j=1:nv
    line1=[num2str(j) ' ' node_name{j} ' ' num2str(xcoord(j)) ' ' ...
        num2str(ycoord(j)) ' ' num2str(p_min(j)) ' ' num2str(p_max(j)) ' ' ...
        num2str(q_min(j)) ' ' num2str(q_max(j)) ' ' num2str(isslack(j)) '\n'];
    fprintf(fid,line1);
end
for j=1:ne
    line1=[num2str(j) ' ' pipe_name{j} ' ' num2str(from_id(j)) ' ' num2str(to_id(j)) ' ' ...
        num2str(diameter(j)) ' ' num2str(pipe_length(j)) ' ' num2str(lambda(j)) ...
        ' ' num2str(disc_seg(j)) '\n'];
    fprintf(fid,line1);
end
for j=1:nc
    line1=[num2str(j) ' ' comp_name{j} ' ' num2str(loc_node(j)) ' ' num2str(to_edge(j)) ' ' ...
        num2str(c_min(j)) ' ' num2str(c_max(j)) ' ' num2str(hp_max(j)) ' ' ...
        num2str(flow_min(j)) ' ' num2str(flow_max(j)) '\n'];
    fprintf(fid,line1);
end
for j=1:ng
    line1=[num2str(j) ' ' gnode_name{j} ' ' num2str(phys_node(j)) '\n'];
    fprintf(fid,line1);
end

fclose('all');