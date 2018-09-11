function [network_out]=gas_model_reconstruct_new(network_in,lmax,disc_opt)
% Anatoly Zlotnik, revised July 2017
%
% includes compressors sorted along link order

% usage: [network]=gas_model_reconstruct_new(network_in,lmax)
%
%   disc_opt=0  = use disc_seg data
%   disc_opt=1  = override and use lmax only
%
% output: network struct
%         network_out.nv = number of vertices
%         network_out.ne = number of edges
%         network_out.nc = number of compressors
%         network_out.ng = number of gnodes
% node data:
%         network_out.node_name = node name strings
%         network_out.xcoord = x coordinate
%         network_out.ycoord = y coordinate
%         network_out.p_min = min pressure at vertices (Pa)
%         network_out.p_max = max pressure at vertices (Pa)
%         network_out.q_min = min injection (negative = consumption) at vertices (kg/s)
%         network_out.q_max = max injection (negative = consumption) at vertices (kg/s)
%         network_out.isslack = 1 -> slack node, 0 -> not slack
%         network_out.nomflow = nominal flow, if 9th column in file
% edge data:
%         network_out.pipe_name = pipe name strings
%         network_out.from_id = start node #
%         network_out.to_id = end node #
%         network_out.diameter = pipe diameter (meters \ inches)
%         network_out.pipe_length = pipe length (km \ miles)
%         network_out.lambda = multiplier
%         network_out.disc_seg = segments for discretization by pipe (0 = lmax)
% compressor data:
%         network_out.comp_name = comp name strings
%         network_out.loc_node = compressor location node #
%         network_out.to_edge = edge # to which flow is boosted 
%         network_out.c_min = min compression ratio
%         network_out.c_max = max compression ratio     
%         network_out.hp_max = max compressor power (hp)
%         network_out.flow_min = min flow through compressor (kg/s)
%         network_out.flow_max = max flow through compressor (kg/s)
% gnode data:
%         network_out.gnode_name = gnode name strings
%         network_out.phys_node = gnode physical location
% reconstruction data:
%         network_out.from_flows=from_flows;  %new edge flow index entering old edges
%         network_out.to_flows=to_flows;      %new edge flow index leaving old edges      
%


nv=network_in.nv;              %number of vertices
ne=network_in.ne;              %number of edges
nc=network_in.nc;              %number of comps
ng=network_in.ng;              %number of gnodes
xcoord=network_in.xcoord;            %lattitude at vertices
ycoord=network_in.ycoord;            %longitude at vertices
p_min=network_in.p_min;        %min pressure at vertices (Pa)
p_max=network_in.p_max;        %max pressure at vertices (Pa)
q_min=network_in.q_min;        %min production (consumption) at vertices (kg/s)
q_max=network_in.q_max;        %max production (consumption) at vertices (kg/s)
isslack=network_in.isslack;    %slack node bool
from_id=network_in.from_id;    %start node #
to_id=network_in.to_id;        %end node #
diameter=network_in.diameter;  %pipe diameter (m)
pipe_length=network_in.pipe_length;      %pipe length (km)
lambda=network_in.lambda;        %multiplier
disc_seg=network_in.disc_seg;   %segments for discretization by pipe (0 = lmax)
loc_node=network_in.loc_node;  %multiplier
to_edge=network_in.to_edge;     %edge # to which flow is boosted 
c_min=network_in.c_min;        %min compression ratio
c_max=network_in.c_max;        %max compression ratio
hp_max=network_in.hp_max;      %max compressor power (hp)
flow_min=network_in.flow_min;      %max compressor power (kg/s)
flow_max=network_in.flow_max;      %max compressor power (kg/s)

%finer discretization of long edges
if(disc_opt==1)
    long_pipes=find(pipe_length>lmax);      %pipes longer than lmax
    seg_div=floor(pipe_length(long_pipes)/lmax)+1;  %how many segments to divide pipes into
elseif(disc_opt==0)
    lpipes=(pipe_length>lmax & disc_seg==0); spipes=(disc_seg>1);
    long_pipes=[find(lpipes);find(spipes)];      %pipes longer than lmax
    %how many segments to divide pipes into
    seg_div=[floor(pipe_length(lpipes)/lmax)+1;disc_seg(spipes)];
end
node_counter=nv;                        %start counting added nodes from nv
edge_counter=ne;                        %start counting added edges from ne
to_edge_new=zeros(nc,1);                %new to edge data for compressors
from_flows=[1:ne]';                 %new edge flow index entering old edges
to_flows=[1:ne]';                   %new edge flow index leaving old edges
for j=1:length(long_pipes)
    newnodes=[node_counter+1:node_counter+seg_div(j)-1]';   %indices of added nodes for long pipe j
    newedges=[edge_counter+1:edge_counter+seg_div(j)]';     %indices of added edges for long pipe j
    from_flows(long_pipes(j))=edge_counter+1-length(long_pipes);
    to_flows(long_pipes(j))=edge_counter+seg_div(j)-length(long_pipes);
    %additional node values
    for k=1:length(newnodes)
        labels{newnodes(k)}='None';
    end
    xcoord(newnodes,:)=(xcoord(from_id(long_pipes(j)))*([seg_div(j)-1:-1:1]/(seg_div(j)))+...
        xcoord(to_id(long_pipes(j)))*([1:seg_div(j)-1]/(seg_div(j))));
    ycoord(newnodes,:)=(ycoord(from_id(long_pipes(j)))*([seg_div(j)-1:-1:1]/(seg_div(j)))+...
        ycoord(to_id(long_pipes(j)))*([1:seg_div(j)-1]/(seg_div(j))));
    p_min1=min([p_min(from_id(long_pipes(j))) p_min(to_id(long_pipes(j)))]);
    p_max1=max([p_max(from_id(long_pipes(j))) p_max(to_id(long_pipes(j)))]);
    p_min(newnodes,:)=p_min1*ones(seg_div(j)-1,1);
    p_max(newnodes,:)=p_max1*ones(seg_div(j)-1,1);
    q_min(newnodes,:)=zeros(seg_div(j)-1,1);
    q_max(newnodes,:)=zeros(seg_div(j)-1,1);
    isslack(newnodes,:)=zeros(seg_div(j)-1,1);
    %additional link values
    diameter(newedges,:)=diameter(long_pipes(j))*ones(seg_div(j),1);
    pipe_length(newedges,:)=pipe_length(long_pipes(j))/seg_div(j)*ones(seg_div(j),1);
    lambda(newedges,:)=lambda(long_pipes(j))*ones(seg_div(j),1);
    from_id(newedges,:)=[from_id(long_pipes(j)); newnodes];
    to_id(newedges,:)=[newnodes;to_id(long_pipes(j))];
    %adjust compressor to_edge if needed
    for k=1:nc
        if(to_edge(k)==long_pipes(j))
            to_edge_new(k)=newedges(1)-length(long_pipes);
        end
    end
    node_counter=node_counter+seg_div(j)-1;     %increase node counter by added nodes
    edge_counter=edge_counter+seg_div(j);       %increase edge counter by added edges
end
%adjust compressor to_edge
for k=1:nc
    if(isempty(intersect(long_pipes,to_edge(k))))
        to_edge_new(k)=to_edge(k)-sum(long_pipes<to_edge(k));
    end
end
%adjust from and to edge of not long_pipes
short_pipes=setdiff([1:ne]',long_pipes);
for k=1:length(short_pipes)
    from_flows(short_pipes(k))=from_flows(short_pipes(k))-sum(long_pipes<short_pipes(k));
    to_flows(short_pipes(k))=to_flows(short_pipes(k))-sum(long_pipes<short_pipes(k));
end
%remove old data for long pipes
diameter(long_pipes,:)=[];
pipe_length(long_pipes,:)=[];
lambda(long_pipes,:)=[];
from_id(long_pipes,:)=[];
to_id(long_pipes,:)=[];

%new network size
nv=node_counter;
ne=edge_counter-length(long_pipes);

%adjacency matrix
Amm=sparse(nv,ne); Amp=sparse(nv,ne);
for i=1:ne
    Amm(from_id(i),i)=-1; Amp(to_id(i),i)=1;
end
Am=Amp+Amm;
% for i=1:nv
%    Ad(i,:)=-(from_id==i)+(to_id==i);
% end
%slack and demand submatrices
As=Am(isslack==1,:);
Ad=Am(isslack==0,:);

%compressor positions
% the array comp_pos has node where compressor is located in first column,
% and the edge that it is compressing into in the second column
comp_pos=[loc_node to_edge_new];
%comp_pos=sortrows(comp_pos,1);

%source nodes
source_nodes=find(q_min>=0 & q_max>0);  %gas source nodes
demand_nodes=find(q_max<=0);            %any node that is not a source
consumer_nodes=find(q_min<0);           %any node that consumes gas
transit_nodes=find(q_max==0 & q_min==0);%any node where no gas is injected or removed from the network
slack_nodes=find(isslack);      %slack nodes
nonslack_nodes=setdiff([1:nv]',slack_nodes); %nodes that are not the slack bus
   
network_out.nv=nv;              %number of vertices
network_out.ne=ne;              %number of edges
network_out.nc=nc;              %number of compressors
network_out.ng=ng;              %number of gnodes
network_out.xcoord=xcoord;      %x coordinate of vertices
network_out.ycoord=ycoord;      %y coordinate of vertices
network_out.p_min=p_min;        %min pressure at vertices (Pa)
network_out.p_max=p_max;        %max pressure at vertices (Pa)
network_out.q_min=q_min;        %min production (consumption) at vertices (kg/s)
network_out.q_max=q_max;        %max production (consumption) at vertices (kg/s)
network_out.isslack=isslack;    %boolean
network_out.from_id=from_id;    %start node #
network_out.to_id=to_id;        %end node #
network_out.diameter=diameter;  %pipe diameter (m)
network_out.pipe_length=pipe_length;      %pipe length (km)
network_out.lambda=lambda;          %multiplier
network_out.loc_node=loc_node;  %multiplier
network_out.to_edge=to_edge_new;    %edge # to which flow is boosted 
network_out.c_min=c_min;        %min compression ratio
network_out.c_max=c_max;        %max compression ratio
network_out.hp_max=hp_max;      %max compressor power (hp)
network_out.flow_min=flow_min;      %min compressor flow (kg/s)
network_out.flow_max=flow_max;      %max compressor flow (kg/s)
network_out.Am=Am;              %adjacency matrix
network_out.As=As;              %slack node adjacency submatrix
network_out.Ad=Ad;              %slack node adjacency submatrix
network_out.Amp=Amp;                %adjacency matrix positive part
network_out.Amm=Amm;                %adjacency matrix negative part
network_out.comp_pos=comp_pos;  %compressor positions (to sub into adjacency matrix)
network_out.source_nodes=source_nodes;  %source nodes
network_out.demand_nodes=demand_nodes;  %demand nodes
network_out.consumer_nodes=consumer_nodes;  %nodes with consumption
network_out.transit_nodes=transit_nodes;    %nodes without consumption
network_out.slack_nodes=slack_nodes;    %slack nodes
network_out.nonslack_nodes=nonslack_nodes;  %non-slack nodes
network_out.phys_node=network_in.phys_node;    %gnode physical node locations
network_out.from_flows=from_flows;  %new edge flow index entering old edges
network_out.to_flows=to_flows;      %new edge flow index leaving old edges                 