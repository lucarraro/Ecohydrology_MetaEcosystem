function [biomassTL, shannonTL] = eval_biomassTL(yy, Amat, bm)
Aplus=Amat; Aplus(Aplus<0)=0;
nNodes=size(yy,2);
biomassTL=zeros(nNodes,100); shannonTL=zeros(nNodes,100); 
for i = 1:nNodes
    y1=yy(1:100,i); % take out D and N?
    z1=y1.*bm; 
    tmp=Aplus.*repmat(y1',100,1);
    M1=tmp./sum(tmp,2);
    M1(isnan(M1))=0;

    M2=eval_TL(M1);
    M3=sum(M2.*z1);% biomass [kg m-3] at each TL
    biomassTL(i,1:length(M3))=M3; 
    % shannon entropy at each TL
    M4=M2.*z1;
    M4=M4./sum(M4);
    LM4=log(M4);
    LM4(LM4==-Inf)=0;
    P=-sum(M4.*LM4);
    shannonTL(i,1:length(P))=P;
end
ind=find(sum(biomassTL)>0); ind=ind(end);
biomassTL=biomassTL(:,1:ind);
shannonTL=shannonTL(:,1:ind);
end