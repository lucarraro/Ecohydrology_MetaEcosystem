function M2=eval_TL(M1)
nSp=size(M1,1);
D = M1>0;
Dtmp = D;
M2 = zeros(nSp,1);

trashbin=[];
noDiet=find(sum(Dtmp,2)==0);
trashbin=[trashbin noDiet];

Dtmp(:,noDiet)=0;
M2(noDiet,1)=1;

while ~isempty(noDiet)
    noDiet=setdiff(find(sum(Dtmp,2)==0), trashbin);
    trashbin=[trashbin; noDiet];
    M3 = [zeros(nSp,1) M2];
    M2 = [M2 zeros(nSp,1)];
    M2(noDiet,:)=M1(noDiet,:)*M3;
    Dtmp(:,noDiet)=0;
end

if (sum(M2(:,end)==0)); M2=M2(:,1:end-1); end

end

