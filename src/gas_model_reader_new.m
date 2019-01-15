function [network]=gas_model_reader_new(mfolder)
% Anatoly Zlotnik, revised August 2017
%
% usage: [network]=gas_model_reader(mod_num)
%
% input: mfolder = name of folder containing network.txt
%
% input file format: space separated file
% first row: [# nodes] [# edges] [# compressors] [#gnodes]
% node rows: [node #] [node name] [xcoord] [ycoord] [pmin] [pmax] [qmin] [qmax] [slack bool]
% edge rows: [edge #] [pipe name] [from node] [to node] [diameter (m)] [length (km)] [friction factor] [disc seg]
% comp rows: [comp #] [comp name] [loc node] [to edge] [cmin] [cmax] [hp max] [min flow] [max flow]
% gnode rows: [gnode #] [gnode name] [phys loc]
% note: numbering begins from 1
%
% output: network struct
%         network.nv = number of vertices
%         network.ne = number of edges
%         network.nc = number of compressors
%         network.ng = number of compressors
% node data:
%         network.node_name = node name strings
%         network.lat = lattitude
%         network.lon = longitude
%         network.p_min = min pressure at vertices (MPa \ psi)
%         network.p_max = max pressure at vertices (MPa \ psi)
%         network.q_min = min injection (negative = consumption) at vertices (kg/sec)
%         network.q_max = max injection (negative = consumption) at vertices (kg/sec)
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
%         network.flow_min = min flow through compressor (kg/s)
%         network.flow_max = max flow through compressor (kg/s)
% gnode data:
%         network.gnode_name = gnode name strings
%         network.phys_node = gnode physical location

fid = fopen([mfolder '\input_network.txt']);
tline = fgetl(fid);
sizes = str2num(tline);

nv=sizes(1);         %number of vertices
ne=sizes(2);         %number of edges
nc=sizes(3);         %number of compressors
ng=sizes(4);         %number of gnodes

%information for each node
xcoord=zeros(nv,1);    %lattitude
ycoord=zeros(nv,1);    %longitude
p_min=zeros(nv,1);  %min pressure
p_max=zeros(nv,1);  %max pressure
q_min=zeros(nv,1);  %min inflow to node (negative for sinks)
q_max=zeros(nv,1);  %max inflow to node (zero for sinks)
isslack=zeros(nv,1);   %1 if slack bus otherwise zero
nomflow=0;   %nominal injection if available
for j=1:nv
   tline = fgetl(fid);
   space_locs=[];
   space_count=0;
   for tcount=1:size(tline,2)
       if(tline(tcount)==' ')
           space_count=space_count+1;
           space_locs(space_count)=tcount;
       end
   end
   node_name{j}=tline(space_locs(1)+1:space_locs(2)-1);
   xcoord(j)=str2num(tline(space_locs(2)+1:space_locs(3)-1));
   ycoord(j)=str2num(tline(space_locs(3)+1:space_locs(4)-1));
   p_min(j)=str2num(tline(space_locs(4)+1:space_locs(5)-1));
   p_max(j)=str2num(tline(space_locs(5)+1:space_locs(6)-1));
   flow_str1=tline(space_locs(6)+1:space_locs(7)-1);
   if(strcmp(flow_str1,'-Infinity')), q_min(j)=-Inf; else
   q_min(j)=str2num(tline(space_locs(6)+1:space_locs(7)-1)); end
   q_max(j)=str2num(tline(space_locs(7)+1:space_locs(8)-1));
   isslack(j)=str2num(tline(space_locs(8)+1:size(tline,2)));
end

%information for each pipe
from_id=zeros(ne,1);    %from node
to_id=zeros(ne,1);      %to node
diameter=zeros(ne,1);   %diameter (m)
pipe_length=zeros(ne,1); %length (km)
lambda=zeros(ne,1);       %friction coefficient 
disc_seg=zeros(ne,1);   %discretization segments (0 = use lmax)
for j=1:ne
   tline = fgetl(fid);
   space_locs=[];
   space_count=0;
   for tcount=1:size(tline,2)
       if(tline(tcount)==' ')
           space_count=space_count+1;
           space_locs(space_count)=tcount;
       end
   end
   pipe_name{j}=tline(space_locs(1)+1:space_locs(2)-1);
   from_id(j)=str2num(tline(space_locs(2)+1:space_locs(3)-1));
   to_id(j)=str2num(tline(space_locs(3)+1:space_locs(4)-1));
   diameter(j)=str2num(tline(space_locs(4)+1:space_locs(5)-1));
   pipe_length(j)=str2num(tline(space_locs(5)+1:space_locs(6)-1));
   lambda(j)=str2num(tline(space_locs(6)+1:space_locs(7)-1));
   disc_seg(j)=str2num(tline(space_locs(7)+1:size(tline,2)));
end

%information for each compressor
%         network.loc_node = compressor location node #
%         network.to_edge = edge # to which flow is boosted 
%         network.c_min = min compression ratio
%         network.c_max = max compression ratio     
%         network.hp_max = max compressor power (hp)
%         network.flow_min = min flow (kg/s)
%         network.flow_max = max flow (kg/s)
%         network.comp_bool = 1 if compressor on edge
loc_node=zeros(nc,1);    %location node for compressor
to_edge=zeros(nc,1);     %to edge
c_min=zeros(nc,1);       %min compression ratio
c_max=zeros(nc,1);       %max compression ratio
hp_max=zeros(nc,1);      %max compressor power (hp)
flow_min=zeros(nc,1);      %max compressor power (hp)
flow_max=zeros(nc,1);      %max compressor power (hp)
comp_bool=zeros(ne,1);   % = 1 if compressor on the edge, 0 otherwise
for j=1:nc
      tline = fgetl(fid);
   space_locs=[];
   space_count=0;
   for tcount=1:size(tline,2)
       if(tline(tcount)==' ')
           space_count=space_count+1;
           space_locs(space_count)=tcount;
       end
   end
comp_name{j}=tline(space_locs(1)+1:space_locs(2)-1);
loc_node(j)=str2num(tline(space_locs(2)+1:space_locs(3)-1));
to_edge(j)=str2num(tline(space_locs(3)+1:space_locs(4)-1));
c_min(j)=str2num(tline(space_locs(4)+1:space_locs(5)-1));
c_max(j)=str2num(tline(space_locs(5)+1:space_locs(6)-1));
hp_max(j)=str2num(tline(space_locs(6)+1:space_locs(7)-1));
flow_min(j)=str2num(tline(space_locs(7)+1:space_locs(8)-1));
flow_max(j)=str2num(tline(space_locs(8)+1:size(tline,2)));
comp_bool(to_edge(j))=1;
end

%information for each gnode
%         network.phys_node = physical node location of gnode
phys_node=zeros(ng,1);    %location node for compressor
for j=1:ng
   tline = fgetl(fid);
   space_locs=[];
   space_count=0;
   for tcount=1:size(tline,2)
       if(tline(tcount)==' ')
           space_count=space_count+1;
           space_locs(space_count)=tcount;
       end
   end
gnode_name{j}=tline(space_locs(1)+1:space_locs(2)-1);
phys_node(j)=str2num(tline(space_locs(2)+1:size(tline,2)));
end

fclose(fid);

comp_pos=[loc_node to_edge];

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
network.pipe_name=pipe_name;      %pipe name strings
network.comp_name=comp_name;      %comp name strings
network.gnode_name=gnode_name;    %gnode name strings
network.xcoord=xcoord;            %lattitude of vertices
network.ycoord=ycoord;            %longitude of vertices
network.p_min=p_min;        %min pressure at vertices (Pa)
network.p_max=p_max;        %max pressure at vertices (Pa)
network.q_min=q_min;        %min injection (negative = consumption) at vertices (kg/sec)
network.q_max=q_max;        %max injection (negative = consumption) at vertices (kg/sec)
network.isslack=isslack;
network.from_id=from_id;    %start node #
network.to_id=to_id;        %end node #
network.diameter=diameter;  %pipe diameter (m)
network.pipe_length=pipe_length;      %pipe length (km)
network.lambda=lambda;      %friction factor
network.disc_seg=disc_seg;  %segments for discretization by pipe (0 = lmax)
network.c_min=c_min;        %min compression ratio
network.c_max=c_max;        %max compression ratio
network.loc_node=loc_node;  %compressor location node #
network.to_edge=to_edge;    %edge # to which flow is boosted 
network.c_min=c_min;        %min compression ratio
network.c_max=c_max;        %max compression ratio     
network.hp_max=hp_max;      %max compressor power (hp)
network.flow_min=flow_min;  %min compressor flow (kg/sec)
network.flow_max=flow_max;  %max compressor flow (kg/sec)
network.nomflow=nomflow;    %nominal injection
network.comp_pos=comp_pos;  %compressor positions (to sub into adjacency matrix)
network.comp_bool=comp_bool;    % = 1 if compressor on edge, = 0 otherwise
network.Am=Am;              %adjacency matrix
network.Amp=Amp;            %adjacency matrix positive part
network.Amm=Amm;            %adjacency matrix negative part
network.slack_nodes=slack_nodes; %slack nodes
network.nonslack_nodes=nonslack_nodes; %nonslack nodes
network.phys_node=phys_node;    %gnode physical locations