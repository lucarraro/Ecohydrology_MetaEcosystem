clear all; close all; clc

OCN_type=144;
scenario='_unifN'; %  ''; '_noLight'; '_unifN'

addpath('utilities')
load('30FW_DF_2_5_8.mat') 
load(['OCN_',num2str(OCN_type),'.mat'])

nSpecies=length(r_list.w1)+2;
nFW=length(fieldnames(r_list));
L(L==0)=100; % fix L=0 at the outlet
nReach=length(downNode);
Q=width.*depth.*velocity;
V=width.*depth.*L;
Wt=W';
if strcmp(scenario,'_noLight')
    light=ones(nReach,1);
else
    light=(1-exp(-depth))./depth.*(1-exp(-width/5)); 
    light=light/mean(light); 
end

kN_vec=[1e-1]; % m-2/day
dispRate_vec=[1e-7]*1e3; % m/s (or 1e-7 s-1 for a length of 1e3 m) 
downBias_vec=[1e1];
vUptake_vec=[1e-5]; % upatake velocity [m/s]

if strcmp(scenario,'_unifN')
    LUscore_mat=ones(nReach,1);
else
    LUscore_mat=1-distToOutlet/2.5e5; 
    LUscore_mat(LUscore_mat<0)=0;
end

detritusRelVel=0.1;
epsilonMineralization=0.1;
epsilonRecycling=1; 
epsilonTerrDetritus=0.1;
bodymass_DN=[4e-5 2e-5];

for ind_FW=1:nFW
    a_mat = A_list.(['w',num2str(ind_FW)]);
    bodymass = bodymass_list.(['w',num2str(ind_FW)]);
    r = r_list.(['w',num2str(ind_FW)]);
    rng(ind_FW); 
    fnam=['r_',num2str(ind_FW),'_',num2str(OCN_type),scenario,'.mat'];    
    if not(isfile(['results/',fnam]))
        tmp=[];
        save(['results/',fnam],'tmp')
        
        kN_mean=kN_vec;
        LUscore=LUscore_mat;
        meanDispRate=dispRate_vec;
        downBias=downBias_vec;
        vUptake=vUptake_vec;
        
        kN=kN_mean*LUscore*max(A)/sum(As.*LUscore)/86400;
       
        dispersalRate=nSpecies*meanDispRate*(bodymass.^0.36)./sum(bodymass.^0.36);
        pD_vec = 0.5*(1 + exp(-downBias*bodymass));
        p_matrix=eval_p_matrix(nReach,pD_vec,downNode,V,depth,width); 
        
        parameters=v2struct(Q,V,L,Wt,kN,kN_mean,As,r,a_mat,vUptake,depth,p_matrix,...
            dispersalRate,light,detritusRelVel,epsilonMineralization,epsilonRecycling, ...
            epsilonTerrDetritus,bodymass,bodymass_DN);
        tic;
        y_mat=zeros(nSpecies*nReach,1); y_mat(:,1)=1e-4*ones(nSpecies*nReach,1);
        ind_D=(nSpecies-1):nSpecies:nReach*nSpecies;
        ind_N=nSpecies:nSpecies:nReach*nSpecies;
        y_mat(ind_D,1)=1000; y_mat(ind_N,1)=100;

        for ind_time=1:100 % in tens of days
            y0=y_mat(:,ind_time);
            tspan=[1+10*(ind_time-1)*86400:100:86400*10*ind_time]; % 10 days
            [t,y] = ODE_ME(parameters,tspan,y0);
            y0=y(end,:)';
            y_mat(:,ind_time+1)=y0;
            fprintf('OCN: %d  -  FW: %d  -  scenario: %s   -   Elapsed time: %.2f s -  Sim time: %d d\n',...
                OCN_type,ind_FW,scenario,toc,10*ind_time)
        end
        y=y_mat(:,end);
        timeSim=1*ind_time; timeElapsed=toc;
        save(['results_newEdge/',fnam],'y','timeSim','timeElapsed',...
            'kN_vec','dispRate_vec','downBias_vec','vUptake_vec','LUscore_mat','bodymass')
    else
        fprintf('FW: %d  -  done!\n',ind_FW)
    end
    disp(' ')
end



