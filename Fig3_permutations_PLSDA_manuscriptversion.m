function datamean = Fig3_permutations_PLSDA_manuscriptversion(filename)

%% Preliminaries
% clc; clear all; close all;
clear;
clc;
% make sure matlab has folder in its path with all of the figure scripts!
% must have bioinformatics plugin installed to run
% addpath('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code');
%% Load in the data
% dataPath='cytokineDataCorrected.xlsx';
dataPathanalyte='cytokineDataCorrected.xlsx';

% dataPath = filename;

dataPath='C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Figure3_Code\cytokineDataCorrected.xlsx';
% dataPath='C:\Users\kgarza6\Box\Project_FlickerNeuroimmune_Team\Data\KristieRawFlickerData\Flicker 701-760_microgliadepeltion_Nov2020\cytokineDataCorrected.xlsx';

% [num txt raw]=xlsread(dataPath,'sheet1','D1:AI1');
[num txt raw]=xlsread(dataPathanalyte,'sheet1','D1:AI1');
analyteLabels=txt;

[num txt raw]=xlsread(dataPath,'sheet1', 'D2:AI64'); 
MFIval=num;

[num txt raw]=xlsread(dataPath, 'sheet1','C2:C64');
sampleLabels=raw;
% 
[num txt raw]=xlsread(dataPath, 'sheet1','B2:B64');
GroupName=txt;

% AMP 6/6/22 added to color correct images and error ellipse
othercolor = ('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\othercolor');
error_ellipse = ('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\error_ellipse');
addpath(othercolor);
addpath(error_ellipse);
%% INPUT:
%name of background wells:
backgroundname = 'Background';

%groupnames:
GroupsforAnalysis=["40Hz";"40HzPLX";"20Hz";"NSPLX";"NS"];
% GroupsforAnalysis=["40Hz";"40HzPLX";"20Hz"];
 

%%
% %Find background wells and average replicates
backgroundrows =[];
counter = 1;
for i = 1:length(sampleLabels)
    if strcmp(backgroundname,sampleLabels(i))
        backgroundrows(counter) = i;
        counter = counter + 1;
    end
end
background = mean(MFIval(backgroundrows,:));

%%
%%Find Analysis Groups by comparing name of groups to group names in excel
%%doc
cntA=1;
indgrouprows=[];
for i=1:length(GroupName);
    for j=1:length(GroupsforAnalysis);
        if strcmp(GroupsforAnalysis(j),GroupName(i));
            indgrouprows(cntA)=i;
             cntA=cntA+1;
   end
    end
end
indgrouprows=indgrouprows';

% only keep analysis groups
 sampleLabels=sampleLabels(indgrouprows);
 MFIval=MFIval(indgrouprows,:);
 GroupName=GroupName(indgrouprows);
 
%reorder
[GroupName,groupsordered]=sort(GroupName);
indgrouprows=indgrouprows(groupsordered);
sampleLabels=sampleLabels(groupsordered);
MFIval=MFIval(groupsordered,:);

% %% Subtract Background
% %activate this if not normalizing by region
%     for i =1:size(MFIval,2)
%        MFIval(:,i)=MFIval(:,i)-background(i); 
%     end
    MFIval(MFIval<0)=0;
 %%
 
 % Groups
numgroups = length(unique(GroupName));
% 
% %total animals
% 
totalnumsamples=(length(GroupName));

% %Remove groups from earlier GroupName table based on whats left in data MFIval
[x,y]=size(MFIval);
sampleLabels((x+1):end,:)=[];
GroupName((x+1):end,:)=[];

[BarLabels,ta,treatment]=unique(GroupName,'stable');
BarLabels=string(BarLabels);


%%
%Samples

%Preallocate matrices for average mean values and standard deviations
avgMFI = zeros(numgroups,length(analyteLabels));
stdMFI = zeros(numgroups,length(analyteLabels));
anovamatrix={};


%Find Sample Groups and Average
counter3=1;
grouprows={};

for i=1:length(BarLabels);
    grouprows(:,counter3)={find(treatment==i)};
    counter3= counter3+1;
end


%Get averages for each group and analyte
for i=1:length(grouprows);
    for j=1:length(analyteLabels);
    avgMFI(i,j)=mean(MFIval(grouprows{:,i},j));
    end
end

%Get standard deviations for each group and analyte
for i=1:length(BarLabels);
    for j=1:length(analyteLabels);
     stdMFI(i,j)=std(MFIval(grouprows{:,i},j));
    end
end

%Construct ANOVA cell array
for i=1:length(BarLabels);
    for j=1:length(analyteLabels);
        anovamatrix{i,j}=(MFIval(grouprows{:,i},j));
    end
end


%% ANOVA
%table to store p values
pvals=zeros(1,length(analyteLabels));
gnames=BarLabels;


if numgroups>2==1
    for i=1:length(analyteLabels);
    [p,tbl,stats]=anova1(MFIval(:,i),GroupName,'off');
    pvals(1,i)=p;
%     if p<=.05==1
%         figure
%         [c,m,h,gnames]=multcompare(stats,'ctype','bonferroni');
%         title(analyteLabels(i));    
%     end
    end
else
    for j=1:length(analyteLabels);
    [T,p]=ttest2(MFIval((treatment==1),j),MFIval((treatment==2),j));
    pvals(1,j)=p;
    end
end

% % Bar plots for cytokines
for i=1:length(analyteLabels); 
   if i==1 || i==10 || i==19 || i==28 || i==37   
        plotCount=0;
        figure; 
    end
    plotCount=plotCount+1;
    subplot(3,3, plotCount) 
    hold on
    bar(avgMFI(:,i))
    error = stdMFI(:,i)/sqrt(length(grouprows));
    errorbar(avgMFI(:,i),error,'.')
    title(analyteLabels(i));
    xlabel(['p-value= ',num2str(pvals(i))])
    set(gca,'XTick',1:length(BarLabels),'XTickLabel',BarLabels,'XTickLabelRotation',45)
    hold off
end

%%Make datasave- data rearranged in order by group for heatmap and plsr
datasave=MFIval;

%% Heatmap
%Zscore data
zData=zscore(datasave);

%Make heatmap
sampleLabels =sampleLabels;

%% turning on this below 1/10/23
figure('units','normalized','position',[.1 .1 .45 .6])

%Pad the data, so that pcolor will work right
dataFig=zeros(size(datasave)+1); dataFig(1:end-1,1:end-1)=zData;
pcolor(dataFig)
%combine sample name and groupname for y labels
strlabels=string(sampleLabels);
ylabels=[strlabels,GroupName];
ylabels=join(ylabels);
%%flip y 
ylabels=ylabels;
% 
set(gca, 'ytick', [1:size(sampleLabels,1)]+0.5);
set(gca, 'YTickLabel', ylabels(:,:));
set(gca, 'xtick', [1:length(analyteLabels)]+0.5);
set(gca, 'XTickLabel', analyteLabels);
set(gca,'position',[0.1 0.25 0.8 0.65]);

th=rotateticklabel(gca, 90);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
% c=othercolor('RdBu10');
% colormap(flipud(c));
colorbar;
caxis([-3 3]);

%% PLSR

pheno=treatment;
region.pls = pls(zData,treatment, 3,'da');

% Store loadings before the rotation
P1Store=region.pls.P(:,1);
P2Store=region.pls.P(:,2);

% Rotate the first two axes to see the seperating LV more clearly
theta=230; % in degrees
T1Temp=region.pls.T(:,1)*cos(pi/180*theta)+...
    region.pls.T(:,2)*sin(pi/180*theta);
T2Temp=region.pls.T(:,2)*cos(pi/180*theta)+...
    -region.pls.T(:,1)*sin(pi/180*theta);

region.pls.T(:,1)=T1Temp;
region.pls.T(:,2)=T2Temp;

% Rotate the loadings to match the scores
P1Temp=region.pls.P(:,1)*cos(pi/180*theta)+...
    region.pls.P(:,2)*sin(pi/180*theta);

P2Temp=region.pls.P(:,2)*cos(pi/180*theta)+...
    -region.pls.P(:,1)*sin(pi/180*theta);

region.pls.P(:,1)=P1Temp;
region.pls.P(:,2)=P2Temp;


%plot scores LV1-LV2
%//////////////////////////////////////////////////////////////////////
%% Turn ON 1/10/23
% 
fig2=figure('units','normalized','position',[.1 .2 [.2 .25]]);
axes1 = axes('Parent',fig2,'LineWidth',2,'XTick',[-0.5 -0.25 0 0.25 0.5],'FontSize',13,'CLim',[0 3]);
colors = {'b', 'r', [1 0.7 0.5], 'k',[0.3 0.5 0.5],[0.2 0.2 0.2]}';


for i = 1:length(BarLabels)
    ind = find(treatment==i);
    plot(region.pls.T(ind,1),region.pls.T(ind,2), 'o', 'MarkerFaceColor', colors{i,1}, 'MarkerEdgeColor', colors{i,1},'MarkerSize',10);
    hold on
end
% 
xlim([-0.75 0.75])
ylim([-0.75 0.75])
box on
xlabel(['Scores on LV1'], 'FontSize', 13);
ylabel(['Scores on LV2'],'FontSize', 13)
legend(BarLabels);
set(gca, 'LineWidth', 1.5,'XTick',[-0.6:0.2:0.6])
title('PLSDA')

%ellipses - use this for permutation averages per group AMP 6/8/22
 % separate LV scores for groups
    GroupLVscores={};
    datamean=[];
    for i=1:length(BarLabels);
            ind= find(treatment==i);
            GroupLVscores{i}=region.pls.T(ind,(1:2));
    end
    
% get mean for each Group LV
    datamean=[];
    for i=1:length(GroupLVscores);
     datamean(i,:)=mean(GroupLVscores{1,i});
    end
    
%plot ellipses by creating covariance matrices
    
for i = 1:length(BarLabels)
    error_ellipse((cov(GroupLVscores{1,i})),datamean(i,:),.68);
%     hold on
end

%% turn on 1/10/23 below
%plot scores LV2-LV3
fig3=figure('units','normalized','position',[.1 .2 [.2 .25]]);
axes1 = axes('Parent',fig3,'LineWidth',2,'XTick',[-0.5 -0.25 0 0.25 0.5],'FontSize',13,'CLim',[0 3]);

for i = 1:length(BarLabels)
    ind = find(treatment==i);
    plot(region.pls.T(ind,2),region.pls.T(ind,3), 'o', 'MarkerFaceColor', colors{i,1}, 'MarkerEdgeColor', colors{i,1});
    hold on
end

xlim([-0.75 0.75]);
ylim([-0.75 0.75]);
box on
xlabel(['Scores on LV2'], 'FontSize', 13);
ylabel(['Scores on LV3'],'FontSize', 13);
legend(BarLabels);
set(gca, 'LineWidth', 1.5,'XTick',[-0.6:0.2:0.6]);
title('PLSDA LV2 LV3');

%plot scores LV1-LV2-LV3
fig4=figure('units','normalized','position',[.1 .2 [.2 .25]]);
axes1 = axes('Parent',fig4,'LineWidth',2,'XTick',[-0.5 -0.25 0 0.25 0.5],'FontSize',13,'CLim',[0 3]);

for i = 1:length(BarLabels)
    ind = find(treatment==i);
    plot3(region.pls.T(ind,1),region.pls.T(ind,2),region.pls.T(ind,3), 'o', 'MarkerFaceColor', colors{i,1}, 'MarkerEdgeColor', colors{i,1});
    hold on
end
% 
xlim([-0.75 0.75]);
ylim([-0.75 0.75]);
box on
xlabel(['Scores on LV1'], 'FontSize', 13);
ylabel(['Scores on LV2'],'FontSize', 13);
zlabel(['Scores on LV3'],'FontSize', 13);
legend(BarLabels)
set(gca, 'LineWidth', 1.5,'XTick',[-0.6:0.2:0.6]);
title('LV1 vs LV2 vs LV3');

% bar plot the loadings on LV1
fig5=figure('units','normalized','position',[.1 .2 [.35 .25]]);

[region.pls.P(:,1) IX1]= sort(region.pls.P(:,1), 'ascend');
bar([1:length(region.pls.P(:,1))], [region.pls.P(:,1)]/max(abs(region.pls.P(:,1))));%, ...
set(gca, 'xtick', [1:length(analyteLabels)]);
set(gca, 'xTickLabel', analyteLabels (IX1));
set(gca,'position',[0.1 0.35 0.8 0.55]);
h=gca;
th=rotateticklabelLV(h, 90);
title('Signals in LV1');
xlim([0.5 size(zData,2)+0.5]); ylim([-1 1]);

% bar plot the loadings on LV2
fig6=figure('units','normalized','position',[.1 .2 [.35 .25]]);

[region.pls.P(:,2) IX2]= sort(region.pls.P(:,2), 'ascend');
bar([1:length(region.pls.P(:,2))], [region.pls.P(:,2)]/max(abs(region.pls.P(:,2))));%, ...
set(gca, 'xtick', [1:length(analyteLabels)]);
set(gca, 'xTickLabel', analyteLabels (IX2));
set(gca,'position',[0.1 0.35 0.8 0.55]);
h=gca;
th=rotateticklabelLV(h, 90);
title('Signals in LV2');
xlim([0.5 size(zData,2)+0.5]); ylim([-1 1]);
% 
figure('units','normalized','position',[.1 .1 .45 .6]);

% bar plot the loadings on LV3
fig7=figure('units','normalized','position',[.1 .2 [.35 .25]]);

[region.pls.P(:,3) IX3]= sort(region.pls.P(:,3), 'ascend');
bar([1:length(region.pls.P(:,3))], [region.pls.P(:,3)]/max(abs(region.pls.P(:,3))));%, ...
set(gca, 'xtick', [1:length(analyteLabels)]);
set(gca, 'xTickLabel', analyteLabels (IX3));
set(gca,'position',[0.1 0.35 0.8 0.55]);
h=gca;
th=rotateticklabelLV(h, 90);
title('Signals in LV3');
xlim([0.5 size(zData,2)+0.5]); ylim([-1 1]);

figure('units','normalized','position',[.1 .1 .45 .6]);

% Replot cytokines in terms of LV1

%Pad the data, so that pcolor will work right
dataFig=zeros(size(datasave)+1); dataFig(1:end-1,1:end-1)=zData(:,IX1);
pcolor(dataFig);
set(gca, 'ytick', [1:size(sampleLabels,1)]+0.5);
set(gca, 'YTickLabel', ylabels);
set(gca, 'xtick', [1:length(analyteLabels)]+0.5);
set(gca, 'XTickLabel', analyteLabels(IX1));
set(gca,'position',[0.1 0.25 0.8 0.65]);

th=rotateticklabel(gca, 90);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
c=othercolor('RdBu10');
colormap(flipud(c));
colorbar;
caxis([-2 2]);

figure('units','normalized','position',[.1 .1 .45 .6]);
% Replot cytokines in terms of LV2

%Pad the data, so that pcolor will work right
dataFig=zeros(size(datasave)+1); dataFig(1:end-1,1:end-1)=zData(:,IX1);
pcolor(dataFig);
set(gca, 'ytick', [1:size(sampleLabels,1)]+0.5);
set(gca, 'YTickLabel', ylabels);
set(gca, 'xtick', [1:length(analyteLabels)]+0.5);
set(gca, 'XTickLabel', analyteLabels(IX2));
set(gca,'position',[0.1 0.25 0.8 0.65]);

th=rotateticklabel(gca, 90);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
c=othercolor('RdBu10');
colormap(flipud(c));
colorbar;
caxis([-2.5 2.5]);

%% LV Plots

%Preallocate matrices for average mean values and standard deviations
avglvscores = zeros(length(BarLabels),3);
stdlvscores = zeros(length(BarLabels),3);
anovamatrixlv1={};

%Find Sample Groups and Average
LVscores=region.pls.T;

%Get averages for each group and analyte
for i=1:length(BarLabels)
    for j=1:3
    avglvscores(i,j)=mean(LVscores(grouprows{:,i},j));
    end
end

%Get standard deviations for each group and analyte
for i=1:length(BarLabels)
    for j=1:3
     stdlvscores(i,j)=std(LVscores(grouprows{:,i},j));
    end
end

%Construct ANOVA cell array
for i=1:length(BarLabels)
    for j=1:3
        anovamatrixlv1{i,j}=(LVscores(grouprows{:,i},j));
    end
end

%% ANOVA
% table to store p values
pvals=zeros(1,3);
gnames=BarLabels;

if numgroups>2==1
    for i=1:3
        [p,tbl,stats]=anova1(LVscores(:,i),GroupName,'off');
        pvals(1,i)=p;
    if p<=.05==1
        figure;
        [c,m,h,gnames]=multcompare(stats,'ctype','bonferroni'); 
        title(i);
    end
    end
else
    for j=1:3
        [T,p]=ttest2(LVscorees((treatment==1),j),LVscores((treatment==2),j));
        pvals(1,j)=p;
    end
end

% % Bar plots for lv % turn On 1/10/23
 figure; 
    title('LV')
for i=1:3 
    subplot(3,1,i);
    hold on
    bar(avglvscores(:,i));
    error = stdlvscores(:,i)/sqrt(length(grouprows));
    errorbar(avglvscores(:,i),error,'.');
    title(i)
    xlabel(['p-value= ',num2str(pvals(i))])
    set(gca,'XTick',1:length(BarLabels),'XTickLabel',BarLabels,'XTickLabelRotation',45)
end

%% Monte Carlo shuffle
% 
% for i=1:1000
%     %Randomly take all but 1 sample data points
%     ind=randperm(size(datasave,1)); 
%     
%     zDataShuffle=zData(ind,:);
%     
%     %Generate PLS Model
%     region.pls = pls(zDataShuffle,treatment,3,'da');
%      
%     tempP1=region.pls.P(:,1);
%     tempP2=region.pls.P(:,2); 
%     
%     tempT1=region.pls.T(:,1);
%     tempT2=region.pls.T(:,2); 
%     
%     %Correct for axis flips in the bootstrap using a scalar product
%     tempP1=sign(sum(tempP1.*P1Store))*[tempP1];
%     tempP2=sign(sum(tempP2.*P2Store))*[tempP2];
% 
%     tempT1=sign(sum(tempP1.*P1Store))*[tempT1];
%     tempT2=sign(sum(tempP2.*P2Store))*[tempT2];
% 
%         
%     %Rotate the loadings to match the scores (same theta as above)
%     rotP1=tempP1*cos(pi/180*theta)+...
%         tempP2*sin(pi/180*theta);
%     rotP2=tempP2*cos(pi/180*theta)+...
%         -tempP1*sin(pi/180*theta);
%     
%     rotT1=tempT1*cos(pi/180*theta)+...
%         tempT2*sin(pi/180*theta);
%     rotT2=tempT2*cos(pi/180*theta)+...
%         -tempT1*sin(pi/180*theta);
%     
%     LV1Shuffle(:,i)=rotP1; 
%     LV2Shuffle(:,i)=rotP2;
%     
%     Scores1Shuffle(:,i)=rotT1; 
%     Scores2Shuffle(:,i)=rotT2;
%      
% end

%Flip LV2 axis to match above
% lv2flip=0;
% if lv2flip
%     LV2_loocv=-LV2_loocv;
% end


%% Monte Carlo sub-sampling, randomly sampling without replacement
% 
% for i=1:1000;
%     %Randomly take all but 1 sample data points
%     ind=randperm(size(datasave,1)); 
%     ind=ind(1:size(datasave,1)-1); %Use this line for LOOCV
%     
%     zDataBoot=zscore(datasave(ind,:));
%     phenoBoot=treatment(ind);
%     
%     
%     %Generate PLS Model
%     region.pls = pls(zDataBoot,(phenoBoot),3,'da');
%      
%     tempP1=region.pls.P(:,1);
%     tempP2=region.pls.P(:,2);  
%     
%     %Correct for axis flips in the bootstrap using a scalar product
%     tempP1=sign(sum(tempP1.*P1Store))*[tempP1];
%     tempP2=sign(sum(tempP2.*P2Store))*[tempP2];
% 
%         
%     %Rotate the loadings to match the scores (same theta as above)
%     rotP1=tempP1*cos(pi/180*theta)+...
%         tempP2*sin(pi/180*theta);
%     rotP2=tempP2*cos(pi/180*theta)+...
%         -tempP1*sin(pi/180*theta);
%     
%     LV1_loocv(:,i)=rotP1; 
%     LV2_loocv(:,i)=rotP2;
%      
% end

%Bar plot mean and std of loadings on LV1 and LV2 for leave-one-out cross
%validation

%Flip LV2 axis to match above
% lv2flip=0;
% if lv2flip
%     LV2_loocv=-LV2_loocv;
% end
% 
% meanLV1=mean(LV1_loocv,2);
% stdLV1=std(LV1_loocv,[],2);
% 
% meanLV2=mean(LV2_loocv,2);
% stdLV2=std(LV2_loocv,[],2);

%Bar plot LV1 and LV2 using the loading ordering (IX1 and IX2) used above.
% fig2=figure('units','normalized','position',[.1 .35 [.3 0.25]]);
% bar([1:length(region.pls.P(:,1))], [meanLV1(IX1)]/max(abs(meanLV1(IX1))));
% hold on
% H=errorbar([1:length(analyteLabels)],[meanLV1(IX1)]/max(abs(meanLV1(IX1))),[stdLV1(IX1)]/max(abs(meanLV1(IX1))), 'k', 'linestyle', 'none', 'linewidth', 1.5);
% set(gca, 'xtick', [1:length(analyteLabels)])
% set(gca, 'xTickLabel', analyteLabels (IX1))
% set(gca,'position',[0.1 0.35 0.8 0.55])
% h=gca;
% th=rotateticklabelLV(h, 90);
% xlim([0.5 size(datasave,2)+0.5]); ylim([-1.25 1.25])
% title('LV1')
% 
% fig3=figure('units','normalized','position',[.1 .35 [.3 0.25]]);
% bar([1:length(region.pls.P(:,1))], [meanLV2(IX2)]/max(abs(meanLV2(IX2))));%, ...
% hold on
% H=errorbar([1:length(analyteLabels)],[meanLV2(IX2)]/max(abs(meanLV2(IX2))),[stdLV2(IX2)]/max(abs(meanLV2(IX2))), 'k', 'linestyle', 'none', 'linewidth', 2);
% set(gca, 'xtick', [1:length(analyteLabels)])
% set(gca, 'xTickLabel', analyteLabels (IX2))
% set(gca,'position',[0.1 0.35 0.8 0.55])
% h=gca;
% th=rotateticklabelLV(h, 90);
% xlim([0.5 size(datasave,2)+0.5]); ylim([-1.25 1.25])
% title('LV2')
% 
% 
% 
% 
% clustergram(zData, 'RowLabels', GroupName, 'ColumnLabels', analyteLabels)

% AMP 6/6/22
disp ('done')
end

