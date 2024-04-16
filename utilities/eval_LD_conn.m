function [linkDensity,connectance] = eval_LD_conn(PA,D)
%EVAL_LD_CONN Summary of this function goes here
%   Detailed explanation goes here
nNodes=size(PA,2);
alpha=sum(PA);
linkDensity=zeros(nNodes,1); connectance=zeros(nNodes,1); 
for i=1:nNodes
    nLinks=sum(sum(D(PA(:,i),PA(:,i))));
    linkDensity(i)=nLinks/alpha(i);
    connectance(i)=nLinks/alpha(i)^2;
end
end

