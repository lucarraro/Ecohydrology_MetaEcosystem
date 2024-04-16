clear all; close all; clc

addpath('utilities')
load('30FW_DF_2_5_8.mat')
bodymass_DN=[4e-5; 2e-5];
exportFigures = false;

%% Read data
% Read OCNs
for OCN_type=[144 480]
    load(['OCN_',num2str(OCN_type),'.mat'])
    L(L==0)=100; % fix L=0 at the outlet
    V=width.*depth.*L;
    nReach=length(downNode);
    switch OCN_type
        case 144
            A144=A;   Z144=Z; L144=L; nReach144=nReach; depth144=depth; width144=width;
            dO144=distToOutlet; As144=As; V144=V; slope144=slope; pathL144=pathLong; pathM144=pathMedian;
            light144=(1-exp(-depth))./depth.*(1-exp(-width/5)); 
        case 480
            A480=A; Z480=Z; L480=L; nReach480=nReach; depth480=depth; width480=width;
            dO480=distToOutlet; As480=As; V480=V; slope480=slope; pathL480=pathLong; pathM480=pathMedian;
            light480=(1-exp(-depth))./depth.*(1-exp(-width/5)); 
    end
end

% initialize variables
biomassTL2=[]; biomassTL3=[]; biomassTL5=[]; biomassTL10=[]; ratioBioTL3=[]; ratioBioTL10=[]; mean_ratioD=[]; shannonH=[]; meanBC=[];  biomassN=[]; biomassD=[];
scenario_vec={'def';'noLight';'unifN'};
str_DN={'N2D8'; 'N5D5'; 'N8D2'};
kolors=[0.66 0.21 0.66; 0.21 0.59 0.21];
for OCN_type=[144 480]
    eval(['nReach=nReach',num2str(OCN_type),';'])
    for NN=[2 5 8]
        DD=10-NN;
        for ind_sc=1:3
            scenario=scenario_vec{ind_sc};
            biomassTL2.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            biomassTL3.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            biomassTL5.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            biomassTL10.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            ratioBioTL3.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            ratioBioTL10.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            shannonTL3.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            shannonTL10.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            totalBiomass.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            biomassD.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            biomassN.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            ratioDN.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
          
            mean_ratioD.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            shannonH.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            meanBC.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            meanBCbal.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            meanBCgra.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            ratioBCbal.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            biomassNF.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
            biomassDF.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)=zeros(nReach,10);
        end
    end
end

% Compute variables
for indFW=1:30
    ind = mod(indFW,10) + 10*(mod(indFW,10)==0);
    for OCN_type=[144 480]
        for ind_sc=1:3
            scenario=scenario_vec{ind_sc};
            if strcmp(scenario,'def')
                scenario_fn='';
            else
                scenario_fn=strcat('_',scenario);
            end
            load(['results\r_',num2str(indFW),'_',num2str(OCN_type),scenario_fn,'.mat'])
            nSpecies=length(r_list.(['w',num2str(indFW)]))+2;
            y_mat=reshape(y,nSpecies,length(y)/nSpecies);
            D=dietMatrix_list.(['w',num2str(indFW)]);
            dFeed=detritusFeeders_list.(['w',num2str(indFW)]);
            nFeed=nutrientFeeders_list.(['w',num2str(indFW)]);
            bm=[bodymass_list.(['w',num2str(indFW)]);bodymass_DN];
            a_mat=A_list.(['w',num2str(indFW)]); a_mat(a_mat<0)=0; % only positive coeffs
            r=r_list.(['w',num2str(indFW)]);
            z_mat=y_mat.*bm;

            [bioTL, shanTL]=eval_biomassTL(y_mat,a_mat,bm);
            ratioBioTL=bioTL./sum(bioTL(:,2:end),2); % biomass of D, N is excluded
            [H, mBC, mBCbal, mBCgra]=eval_alphabeta(y_mat);

            NN=length(nFeed); DD=length(dFeed);
            biomassTL2.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=bioTL(:,2);
            biomassTL3.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=bioTL(:,3);
            biomassTL5.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=bioTL(:,5);
            biomassTL10.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=bioTL(:,10);
            ratioBioTL3.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=ratioBioTL(:,3);
            ratioBioTL10.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=ratioBioTL(:,10);
            shannonTL3.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=shanTL(:,3);
            shannonTL10.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=shanTL(:,10);
            totalBiomass.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=sum(z_mat(1:98,:));
            biomassD.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=z_mat(99,:);
            biomassN.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=z_mat(100,:);
            ratioDN.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=z_mat(99,:)./z_mat(100,:);
     
            mean_ratioD.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=...
                eval_mean_ratioD(D, nSpecies, dFeed, nFeed, a_mat, y_mat, bm);
            shannonH.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=H;
            meanBC.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=mBC;
            meanBCbal.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=mBCbal;
            meanBCgra.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=mBCgra;
            ratioBCbal.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=mBCbal./mBC;
            biomassNF.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=sum(z_mat(nFeed,:),1);
            biomassDF.(['o',num2str(OCN_type)]).(['N',num2str(NN),'D',num2str(DD)]).(scenario)(:,ind)=sum(z_mat(dFeed,:),1);

        end
    end
end


%% Draw main figures

% Fig. 2
f=figure('units','centimeters','position',[0 0 17 18]); 
T=tiledlayout(4,21,'TileSpacing','none','Padding','tight')
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassDF.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassDF.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassDF.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassDF.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 3.3e-5],'ytick',[0:1e-5:3e-5],'xticklabel',[]); box off; title(scenario); yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of D feeders [kg m-3]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(biomassDF.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(biomassDF.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 3.3e-5]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassNF.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassNF.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassNF.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassNF.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'xticklabel',[],'ylim',[0 4.95e-5],'ytick',[0:1.5e-5:4.5e-5]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of N feeders [kg m-3]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(biomassNF.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(biomassNF.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 4.955e-5]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassTL3.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassTL3.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassTL3.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassTL3.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 1.65e-5],'xticklabel',[]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of 3rd TL [kg m-3]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(biomassTL3.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(biomassTL3.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 1.65e-5]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassTL10.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassTL10.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassTL10.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassTL10.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 4.4e-11]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of 10th TL [kg m-3]'); else set(gca,'yticklabel',[]); end 
    if ind_sc==3; lgd = legend('EN-LP','SN-LP','EN-SP','SN-SP'); lgd.Location = 'northwest'; end
    xlabel('Drainage area [m^2]')
    nexttile(T,[1 2])
    y1=(mean(biomassTL10.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(biomassTL10.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 4.4e-11]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
if exportFigures; exportgraphics(f,'Fig2_draft.pdf','ContentType','vector'); end

% Fig. 3
f=figure('units','centimeters','position',[0 0 17 9]); 
T=tiledlayout(2,21,'TileSpacing','none','Padding','tight')
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,ratioBioTL3.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,ratioBioTL3.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,ratioBioTL3.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,ratioBioTL3.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0.1 0.21],'xticklabel',[]); box off; title(scenario);  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of 3rd TL [-]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(ratioBioTL3.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(ratioBioTL3.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0.1 0.21]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,ratioBioTL10.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,ratioBioTL10.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,ratioBioTL10.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,ratioBioTL10.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 8.8e-7]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Ratio of biomass of 10th TL [-]'); else set(gca,'yticklabel',[]); end 
    if ind_sc==3; lgd = legend('EN-LP','SN-LP','EN-SP','SN-SP'); lgd.Location = 'northeast'; end
    xlabel('Drainage area [m^2]')
    nexttile(T,[1 2])
    y1=(mean(ratioBioTL10.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(ratioBioTL10.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0  8.8e-7]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
if exportFigures; saveas(f,'Fig3_draft.png'); exportgraphics(f,'Fig3_draft.pdf','ContentType','vector'); end

% Fig. 4
f=figure('units','centimeters','position',[0 0 17 18]); %kolors=get(gca,'colororder');
T=tiledlayout(4,21,'TileSpacing','none','Padding','tight');
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,mean_ratioD.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,mean_ratioD.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,mean_ratioD.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,mean_ratioD.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    plot([4e6 2.5e9],[0.5 0.5],'--','Color',[0.4 0.4 0.4])
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 1.1],'ytick',[0:0.25:1],'xticklabel',[]); box off; yl=get(gca,'ylim'); title(scenario);
    if ind_sc==1; ylabel('Ratio of energy flux from D [-]'); else set(gca,'yticklabel',[]); end %
    nexttile(T,[1 2])
    y1=(mean(mean_ratioD.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(mean_ratioD.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 1.1]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,shannonH.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,shannonH.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,shannonH.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,shannonH.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[1.5 2.38],'ytick',[1.5:0.2:2.38],'xticklabel',[]); box off; yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Shannon entropy'); else set(gca,'yticklabel',[]); end
    nexttile(T,[1 2])
    y1=(mean(shannonH.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(shannonH.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[1.5 2.38]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,meanBC.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,meanBC.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,meanBC.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,meanBC.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0.15 0.645],'ytick',[0.15:0.15:0.6],'xticklabel',[]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Mean Bray-Curtis dissimilarity'); else set(gca,'yticklabel',[]); end    
    nexttile(T,[1 2])
    y1=(mean(meanBC.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(meanBC.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0.15 0.645]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,ratioBCbal.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,ratioBCbal.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,ratioBCbal.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,ratioBCbal.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 0.88],'ytick',[0:0.2:0.88]); box off; yl=get(gca,'ylim');
    xlabel('Drainage area [m^2]')
    if ind_sc==1; ylabel('% turnover'); else set(gca,'yticklabel',[]); end
    if ind_sc==2; lgd = legend('EN-LP','SN-LP','EN-SP','SN-SP'); lgd.Location = 'northeast'; end
    nexttile(T,[1 2])
    y1=(mean(ratioBCbal.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(ratioBCbal.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 0.66]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
if exportFigures; saveas(f,'Fig4_draft.png'); exportgraphics(f,'Fig4_draft.pdf','ContentType','vector'); end


%% Supplementary Figures
%. Fig. S1bc
load('utilities/zoomin144.mat')
f=figure('units','centimeters','position',[0 0 17 9]); 
T=tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
nexttile;
plot_custom(-dO144,biomassDF.o144.N5D5.def,pp.n4,kolors(1,:),'-.',0,1);
for ind_p=1:11
    if ind_p~=4
plot_custom(-dO144,biomassDF.o144.N5D5.def,pp.(['n',num2str(ind_p)]),kolors(1,:),'--.');
    end
end
set(gca,'tickdir','out','ylim',[0 4.5e-5], ...
    'ytick',[0:1.5e-5:4.5e-5],'xtick',[-2.4e5 -2.3e5 -2.2e5]); box off
ylabel('YDF [kg m-3]'); xlabel('Distance to outlet [m]')
nexttile;
semilogx(A144(pp.n4), biomassDF.o144.N5D5.def(pp.n4,:), '--.','Color',kolors(1,:))
hold on; plot_custom(A144, biomassDF.o144.N5D5.def, pp.n4,kolors(1,:),'-.')
set(gca,'tickdir','out','xlim',[5e6 2.2e8],'ylim',[0 4.5e-5],'ytick',[0:1.5e-5:4.5e-5]); box off
ylabel('YDF [kg m-3]'); xlabel('Drainage Area [m^2]')
if exportFigures; exportgraphics(f,'FigS1bc_draft.pdf','ContentType','vector'); end

% Fig. S2
f=figure('units','centimeters','position',[0 0 17 9]); 
T=tiledlayout(2,21,'TileSpacing','none','Padding','tight')
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassD.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassD.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassD.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassD.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 0.88e-3],'ytick',[0:2e-4:0.8e-3],'xticklabel',[]); box off; title(scenario); yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of D [kg m-3]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(biomassD.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(biomassD.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 0.88e-3]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassN.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassN.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassN.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassN.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 2.2e-3],'ytick',[0:5e-4:2e-3]); box off;  yl=get(gca,'ylim');
        xlabel('Drainage area [m^2]')
    if ind_sc==1; ylabel('Biomass of N [kg m-3]'); else set(gca,'yticklabel',[]); end 
    if ind_sc==3; lgd = legend('EN-LP','SN-LP','EN-SP','SN-SP'); lgd.Location = 'northeast'; end
    nexttile(T,[1 2])
    y1=(mean(biomassN.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(biomassN.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 2.2e-3]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
if exportFigures; exportgraphics(f,'FigS2_draft.pdf','ContentType','vector'); end


%. Fig. S3
bins=[6:0.25:10];
Vfrac144=zeros(length(bins)-1,1); Vfrac480=zeros(length(bins)-1,1);
for bb =1:length(bins)-1
    Vfrac144(bb)=sum(V144(log10(A144)>=bins(bb) & log10(A144)<bins(bb+1)))/sum(V144);
    Vfrac480(bb)=sum(V480(log10(A480)>=bins(bb) & log10(A480)<bins(bb+1)))/sum(V480);
end
f=figure('units','centimeters','position',[0 0 16.7 8]); 
T=tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile;
histogram(log10(A144),[6:0.25:10],'FaceColor',kolors(1,:),'Normalization','probability'); hold on; box off; 
histogram(log10(A480),[6:0.25:10],'FaceColor',kolors(2,:),'Normalization','probability')
set(gca,'tickdir','out'); xlabel('Drainage area [m^2]'); ylabel('Fraction of reaches [-]')
nexttile;
bar(6.125:0.25:9.875,Vfrac144,1,'FaceColor',kolors(1,:),'FaceAlpha',0.5); hold on; box off; 
bar(6.125:0.25:9.875,Vfrac480,1,'FaceColor',kolors(2,:),'FaceAlpha',0.5); 
set(gca,'tickdir','out','xtick',[6:10],'ytick',[0:0.1:0.4]); xlabel('Drainage area [m^2]'); ylabel('Fraction of total water volume [-]')
if exportFigures; exportgraphics(f,'FigS3_draft.pdf','ContentType','vector'); end


% Fig. S4
f=figure('units','centimeters','position',[0 0 17 18]); 
T=tiledlayout(4,21,'TileSpacing','none','Padding','tight')
for ind_sc=1:3
    scenario=str_DN{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassDF.o144.(scenario).def,pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassDF.o480.(scenario).def,pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassDF.o144.(scenario).def,pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassDF.o480.(scenario).def,pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 4.4e-5],'ytick',[0:1e-5:4e-5],'xticklabel',[]); box off; title(scenario); yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of D feeders [kg m-3]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(biomassDF.o144.(scenario).def.*V144)/mean(V144))'; y2=(mean(biomassDF.o480.(scenario).def.*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 4.4e-5]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=str_DN{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassNF.o144.(scenario).def,pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassNF.o480.(scenario).def,pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassNF.o144.(scenario).def,pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassNF.o480.(scenario).def,pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'xticklabel',[],'ylim',[0 6.6e-5]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of N feeders [kg m-3]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(biomassNF.o144.(scenario).def.*V144)/mean(V144))'; y2=(mean(biomassNF.o480.(scenario).def.*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 6.6e-5]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=str_DN{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassTL3.o144.(scenario).def,pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassTL3.o480.(scenario).def,pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassTL3.o144.(scenario).def,pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassTL3.o480.(scenario).def,pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 2.2e-5],'xticklabel',[]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of 3rd TL [kg m-3]'); else set(gca,'yticklabel',[]); end 
    nexttile(T,[1 2])
    y1=(mean(biomassTL3.o144.(scenario).def.*V144)/mean(V144))'; y2=(mean(biomassTL3.o480.(scenario).def.*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 2.2e-5]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=str_DN{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,biomassTL10.o144.(scenario).def,pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,biomassTL10.o480.(scenario).def,pathL480,kolors(2,:),'-.')
    plot_custom(A144,biomassTL10.o144.(scenario).def,pathM144,kolors(1,:),'--.'); plot_custom(A480,biomassTL10.o480.(scenario).def,pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 4.4e-11]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Biomass of 10th TL [kg m-3]'); else set(gca,'yticklabel',[]); end 
    if ind_sc==1; lgd = legend('EN-LP','SN-LP','EN-SP','SN-SP'); lgd.Location = 'northeast'; end
    xlabel('Drainage area [m^2]')
    nexttile(T,[1 2])
    y1=(mean(biomassTL10.o144.(scenario).def.*V144)/mean(V144))'; y2=(mean(biomassTL10.o480.(scenario).def.*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 4.4e-11]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
if exportFigures; saveas(f,'FigS4_draft.png'); exportgraphics(f,'FigS4_draft.pdf','ContentType','vector'); end

% Fig. S5
f=figure('units','centimeters','position',[0 0 17 13.5]); kolors=[0.66 0.21 0.66; 0.21 0.59 0.21];%kolors=get(gca,'colororder');
T=tiledlayout(3,21,'TileSpacing','none','Padding','tight');
for ind_sc=1:3
    scenario=str_DN{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,mean_ratioD.o144.(scenario).def,pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,mean_ratioD.o480.(scenario).def,pathL480,kolors(2,:),'-.')
    plot_custom(A144,mean_ratioD.o144.(scenario).def,pathM144,kolors(1,:),'--.'); plot_custom(A480,mean_ratioD.o480.(scenario).def,pathM480,kolors(2,:),'--.')
    plot([4e6 2.5e9],[0.5 0.5],'--','Color',[0.4 0.4 0.4])
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 1.1],'ytick',[0:0.25:1],'xticklabel',[]); box off; yl=get(gca,'ylim'); title(scenario);
    if ind_sc==1; ylabel('Ratio of energy flux from D [-]'); else set(gca,'yticklabel',[]); end %
    nexttile(T,[1 2])
    y1=(mean(mean_ratioD.o144.(scenario).def.*V144)/mean(V144))'; y2=(mean(mean_ratioD.o480.(scenario).def.*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 1.1]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=str_DN{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,shannonH.o144.(scenario).def,pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,shannonH.o480.(scenario).def,pathL480,kolors(2,:),'-.')
    plot_custom(A144,shannonH.o144.(scenario).def,pathM144,kolors(1,:),'--.'); plot_custom(A480,shannonH.o480.(scenario).def,pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[1.2 2.52],'ytick',[1.2:0.2:2.52],'xticklabel',[]); box off; yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Shannon entropy'); else set(gca,'yticklabel',[]); end
    nexttile(T,[1 2])
    y1=(mean(shannonH.o144.(scenario).def.*V144)/mean(V144))'; y2=(mean(shannonH.o480.(scenario).def.*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[1.2 2.52]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=str_DN{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,meanBC.o144.(scenario).def,pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,meanBC.o480.(scenario).def,pathL480,kolors(2,:),'-.')
    plot_custom(A144,meanBC.o144.(scenario).def,pathM144,kolors(1,:),'--.'); plot_custom(A480,meanBC.o480.(scenario).def,pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0.0 0.88],'ytick',[0.0:0.2:0.8],'xticklabel',[]); box off;  yl=get(gca,'ylim');
    if ind_sc==1; ylabel('Mean Bray-Curtis dissimilarity'); else set(gca,'yticklabel',[]); end    
    nexttile(T,[1 2])
    y1=(mean(meanBC.o144.(scenario).def.*V144)/mean(V144))'; y2=(mean(meanBC.o480.(scenario).def.*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0.0 0.88]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
if exportFigures; saveas(f,'FigS5_draft.png'); exportgraphics(f,'FigS5_draft.pdf','ContentType','vector'); end

% Fig. S6
f=figure('units','centimeters','position',[0 0 17 9]); 
T=tiledlayout(2,21,'TileSpacing','none','Padding','tight')
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,meanBCgra.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,meanBCgra.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,meanBCgra.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,meanBCgra.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 0.66],'ytick',[0:0.2:0.66],'xticklabel',[]); box off; yl=get(gca,'ylim');
    if ind_sc==1; ylabel('BC nestedness'); else set(gca,'yticklabel',[]); end
    nexttile(T,[1 2])
    y1=(mean(meanBCgra.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(meanBCgra.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 0.66]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
for ind_sc=1:3
    scenario=scenario_vec{ind_sc}; nexttile(T,[1 5])
    plot_custom(A144,meanBCbal.o144.N5D5.(scenario),pathL144,kolors(1,:),'-.',0,1); plot_custom(A480,meanBCbal.o480.N5D5.(scenario),pathL480,kolors(2,:),'-.')
    plot_custom(A144,meanBCbal.o144.N5D5.(scenario),pathM144,kolors(1,:),'--.'); plot_custom(A480,meanBCbal.o480.N5D5.(scenario),pathM480,kolors(2,:),'--.')
    set(gca,'tickdir','out','xlim',[4e6 2.5e9],'xtick',[1e7 1e8 1e9],'ylim',[0 0.66],'ytick',[0:0.2:0.66]); box off; yl=get(gca,'ylim');
        xlabel('Drainage area [m^2]')
    if ind_sc==1; ylabel('BC turnover'); else set(gca,'yticklabel',[]); end
    if ind_sc==2; lgd = legend('EN-LP','SN-LP','EN-SP','SN-SP'); lgd.Location = 'northeast'; end
    nexttile(T,[1 2])
    y1=(mean(meanBCbal.o144.N5D5.(scenario).*V144)/mean(V144))'; y2=(mean(meanBCbal.o480.N5D5.(scenario).*V480)/mean(V480))';
    boxplot([y1 y2],'BoxStyle','outline','ColorGroup',[1 2],'Colors',kolors,'Widths',0.5,'Symbol','o'); box off;
    hold on; for i=1:10; plot([1 2],[y1(i) y2(i)],'Color',[0.8 0.8 0.8],'Marker','.'); end
    set(gca,'tickdir','out','ylim',[0 0.66]); axis off
    [h,p]=ttest(y1,y2); if p<0.001; title('***','FontSize',7); elseif p<0.01; title('**','FontSize',7); elseif p<0.05; title('*','FontSize',7); else title('-','FontSize',7); end
end
if exportFigures; exportgraphics(f,'FigS6_draft.pdf','ContentType','vector'); end













