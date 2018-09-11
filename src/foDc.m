function [x,D]=foDc(N)

%x=-([0:N]'-N/2)/(N/2);
x=-([0:N+1]'-(N+1)/2)/((N+1)/2);
D=eye(N+1)-diag(ones(N,1),-1); 
D(1,N+1)=-1;
D=-D*(N/2);

