function []=gas_model_plotter_new(network_in)
%Anatoly Zlotnik, updated September 2017

% usage: [network]=gas_model_plotter(network)
%
% input: network
%
%

nv=network_in.nv;              %number of vertices
ne=network_in.ne;              %number of edges
nc=network_in.nc;              %number of edges
lat=network_in.ycoord;            %lattitude at vertices
lon=network_in.xcoord;            %longitude at vertices
p_min=network_in.p_min;        %min pressure at vertices (psi)
p_max=network_in.p_max;        %max pressure at vertices (psi)
q_min=network_in.q_min;        %min production (consumption) at vertices (cubic ft/day)
q_max=network_in.q_max;        %max production (consumption) at vertices (cubic ft/day)
isslack=network_in.isslack;    %slack node bool
from_id=network_in.from_id;    %start node #
to_id=network_in.to_id;        %end node #
diameter=network_in.diameter;  %pipe diameter (inches)
pipe_length=network_in.pipe_length;      %pipe length (miles)
lambda=network_in.lambda;        %multiplier
loc_node=network_in.loc_node;  %multiplier
to_edge=network_in.to_edge;     %edge # to which flow is boosted 
c_min=network_in.c_min;        %min compression ratio
c_max=network_in.c_max;        %max compression ratio
% hp_max=network_in.hp_max;      %max compressor power (hp)

nodes=[lon lat];
starts=[lon(from_id) lat(from_id)];
ends=[lon(to_id) lat(to_id)];

%transport_nodes=find(q_min==q_max);
%consumer_nodes=find(q_min<0 & q_max<=0);
%source_nodes=find(q_min>=0 & q_max>0);
compressor_links=to_edge;
other_links=setdiff([1:ne],to_edge);
pressure_nodes=find(isslack);
flow_nodes=find(~isslack);

clf
lp1=nodes(pressure_nodes(1),1); lp2=nodes(pressure_nodes(1),2);
plot(lp1,lp2,'b.'), hold on, plot(lp1,lp2,'r.'), %plot(lp1,lp2,'.','Color',[0 .5 0]), %, 
plot(lp1,lp2,'k','LineWidth',3), %plot(lp1,lp2,'.r')
quiver(lp1,lp2,.0001,.0001,'Color','r','MaxHeadSize',5,'LineWidth',3), %plot(lp1,lp2,'r','LineWidth',3), 
legend('flow nodes','pressure nodes','pipes','compressors','Location','NorthWest')
msize=23;
for j=1:nv
    text(nodes(j,1),nodes(j,2)+.03,num2str(j),'Color','b') 
end
plot([starts(other_links,1)';ends(other_links,1)'],...
    [starts(other_links,2)';ends(other_links,2)'],'k','LineWidth',3),
plot([starts(compressor_links,1)';ends(compressor_links,1)'],...
    [starts(compressor_links,2)';ends(compressor_links,2)'],'r','LineWidth',5),
for j=1:nc
    quiver(starts(compressor_links(j),1),starts(compressor_links(j),2),...
        ends(compressor_links(j),1)-starts(compressor_links(j),1),...
        ends(compressor_links(j),2)-starts(compressor_links(j),2),...
        'Color','r','MaxHeadSize',5,'Linewidth',5)
end
 plot(nodes(flow_nodes,1),nodes(flow_nodes,2),'b.','MarkerSize',msize), hold on,
%plot(nodes(consumer_nodes,1),nodes(consumer_nodes,2),'.','Color',[0 .5 0],'MarkerSize',msize)
plot(nodes(pressure_nodes,1),nodes(pressure_nodes,2),'.r','MarkerSize',msize)
for j=1:ne-nc
    xtext=(nodes(from_id(j),1)+nodes(to_id(j),1))/2;
    ytext=(nodes(from_id(j),2)+nodes(to_id(j),2))/2+.03;
    text(xtext,ytext,num2str(j),'Color','k') %,'FontSize',12) 
end
for j=1:nc
    xtext=(nodes(from_id(j+ne-nc),1)+nodes(to_id(j+ne-nc),1))/2;
    ytext=(nodes(from_id(j+ne-nc),2)+nodes(to_id(j+ne-nc),2))/2+.03;
    text(xtext,ytext,num2str(j),'Color','r') %,'FontSize',12) 
end
 hold off