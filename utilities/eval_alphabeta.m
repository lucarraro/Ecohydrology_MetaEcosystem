function [H, meanBC, meanBCbal, meanBCgra, BC] = eval_alphabeta(y_mat)

y_mat=y_mat(1:98,:); % delete D, N

nSp=size(y_mat,1); %98
nReach=size(y_mat,2);
p=y_mat./repmat(sum(y_mat),size(y_mat,1),1); % proportion of individuals

H=zeros(1,size(p,2));
for i=1:size(p,2)
    H(1,i)=-sum(p(:,i).*log(p(:,i))); % Shannon entropy
end

BC=zeros(nReach);  BCbal=BC; BCgra=BC;
for i=1:nReach
    for j = 1:nReach
        BC(i,j)=1-2*sum(min(y_mat(:,i),y_mat(:,j)))/(sum(y_mat(:,i)+y_mat(:,j)));
        % decomposition following Baselga 2013  https://doi.org/10.1111/2041-210X.12029
        A=sum(min(y_mat(:,i),y_mat(:,j)));
        B=sum(y_mat(:,i) - min(y_mat(:,i),y_mat(:,j)));
        C=sum(y_mat(:,j) - min(y_mat(:,i),y_mat(:,j)));
        BCbal(i,j)=min(B,C)/(A+min(B,C));
        BCgra(i,j)=abs(B-C)/(2*A+B+C)*A/(A+min(B,C));
    end
end
meanBC=mean(BC)*(nReach)/(nReach-1); % the factor corrects for the 0s in the diagonal (it should be the mean across nReach-1 values)
meanBCbal=mean(BCbal)*(nReach)/(nReach-1);
meanBCgra=mean(BCgra)*(nReach)/(nReach-1);