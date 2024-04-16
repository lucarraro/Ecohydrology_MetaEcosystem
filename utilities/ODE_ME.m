function [t,y] = ODE_ME(parameters,tspan,y0)

[Q,V,leng,Wt,kN,kN_mean,As,r,a_mat,vUptake,depth,p_matrix,dispersalRate,...
    light,detritusRelVel,epsilonMineralization,epsilonRecycling,epsilonTerrDetritus,...
    bodymass,bodymass_DN]=v2struct(parameters);
nSpecies=size(a_mat,1); nReach=length(Q);

% known term
phi_mat=zeros(nSpecies,nReach);
phi_mat(nSpecies-1,:)=(epsilonTerrDetritus*kN_mean/86400*As./V)';
phi_mat(nSpecies,:)=(kN.*As./V)'; % alternative: (kN_mean*As./V)';
PHI=phi_mat(:);

% loss term
L_mat=zeros(nSpecies,nReach);
L_mat(1:nSpecies-2,:)=repmat(r,1,nReach)-dispersalRate*(1./leng');
L_mat(nSpecies-1,:)=-vUptake./depth'-detritusRelVel*(Q./V)';
L_mat(nSpecies,:)=-vUptake./depth'-(Q./V)';
L=L_mat(:);

% interaction matrix
Amat=sparse(nSpecies*nReach,nSpecies*nReach);
for i=1:nReach
    a_mat2=a_mat;
    ind=(i-1)*nSpecies+1:i*nSpecies;
    a_mat2(end,:) = a_mat(end,:)*light(i); 
    a_mat2(:,end) = a_mat(:,end)*light(i); 
    Amat(ind,ind)=a_mat2;
end

% dispersal matrix
Dmat=sparse(nSpecies*nReach,nSpecies*nReach);
for sp=1:nSpecies-2
    ind_sp=sp:nSpecies:(nSpecies*nReach);
    Dmat(ind_sp,ind_sp)=p_matrix(:,:,sp)*dispersalRate(sp)./repmat(leng',nReach,1); % divide by leng of the origin node to get the right dispRate
end
for i=1:nReach
    Dmat(i*nSpecies-1, (i-1)*nSpecies+(1:nSpecies-2))=-epsilonRecycling*(r.*bodymass)/bodymass_DN(1)';
    Dmat(i*nSpecies, i*nSpecies-1)=epsilonMineralization*vUptake/depth(i)*bodymass_DN(1)/bodymass_DN(2); 
end
ind_D=(nSpecies-1):nSpecies:nReach*nSpecies;
ind_N=nSpecies:nSpecies:nReach*nSpecies;
TransportMatrix=Wt.*repmat(Q',nReach,1); % only upstream transport (don't multiply by 86400 as a's, b's are already in s^-1)
TransportMatrix=sparse(TransportMatrix);
Dmat(ind_D,ind_D)  =  detritusRelVel*1./V.*TransportMatrix;
Dmat(ind_N,ind_N)  =  1./V.*TransportMatrix;
        
   
[t,y]=ode23(@(t,y)eqs(t,y),tspan,y0); %ode23 
%'NonNegative',1:(N_reach*N_var),  ) ,,odeset('RelTol',1e-6,'AbsTol',1e-12) odeset('NonNegative',1:(N_reach*6) 

function dy=eqs(t,y)
    % y is the CONCENTRATION in the reach!! hence the division by V
    dy=zeros(nSpecies*nReach,1);
    
    dy = PHI + L.*y + Dmat*y + (Amat*y).*y;  
    

end

end

