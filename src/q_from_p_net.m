function [q]=q_from_p_net(p,comps,par)
Dk=par.Dk; lamk=par.lamk; ML=par.ML; Ma=par.Ma;
cn=par.n.comp_pos(:,1); cl=par.n.comp_pos(:,2);
spos=par.spos; dpos=par.dpos;
dnodes=par.n.demand_nodes;
snodes=par.n.source_nodes;

   Amj=par.n.Ad; %generate incidence matrix
   for ccc=1:length(dpos)   %add compression values at time t for demand incidence submatrix
       cc=dpos(ccc);
       Amj(cn(cc),cl(cc))=-comps(cc);   
   end
   for ccc=1:length(spos)   %add compression values at time t for slack incidence submatrix
       cc=spos(ccc);
       Amj(cn(cc),cl(cc))=-comps(cc);   
   end
   Adj=Amj(dnodes,:);       %this is demand submatrix
   Asj=Amj(snodes,:);       %this is slack submatrix
   ps=par.n.p_min(snodes)/par.psc;      %this is slack bus density

z=Ma*(Asj'*ps+Adj'*p).*(abs(Asj')*ps+abs(Adj')*p)./(lamk./Dk*par.R);   %flux dynamics
q=sqrt(abs(z)).*sign(z);
