function [mean_ratioD, PP] = eval_mean_ratioD(D, nSpecies, dFeed, nFeed, a_mat, y_mat, bm)

nReach=size(y_mat,2);
Dtmp=D;
trashbin=[];
detritusFlux=zeros(nSpecies,nReach);
nutrientFlux=zeros(nSpecies,nReach);

noDiet=setdiff(find(sum(Dtmp,1)==0), trashbin);
for i = noDiet
    if ismember(i,dFeed)
        detritusFlux(i,:)=a_mat(i,nSpecies-1)*y_mat(nSpecies-1,:)*bm(i);
    end
    if ismember(i,nFeed)
        nutrientFlux(i,:)=a_mat(i,nSpecies)*y_mat(nSpecies,:)*bm(i);
    end
end

while ~isempty(noDiet)
    trashbin=[trashbin noDiet];
    Dtmp(noDiet,:)=0;
    noDiet=setdiff(find(sum(Dtmp,1)==0), trashbin);
    ratioD=detritusFlux./(detritusFlux+nutrientFlux);
    ratioD(nSpecies-1,:)=1; ratioD(nSpecies,:)=0;
    for i = noDiet
        detritusFlux(i,:)=sum(bm(i)*a_mat(i,:)'.*y_mat.*ratioD,1,'omitnan');
        nutrientFlux(i,:)=sum(bm(i)*a_mat(i,:)'.*y_mat.*(1-ratioD),1,'omitnan');
    end
end
ratioD(y_mat==0)=NaN; % ratioD is calculated also for absent species!
mean_ratioD=mean(ratioD,1,'omitnan'); % percentage of energy derived from detritus, averaged across living species (even though the small ones are the most abundant?)

PP=sum(nutrientFlux(nFeed,:).*y_mat(nFeed,:)); % primary production

end