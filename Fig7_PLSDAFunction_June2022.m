%% Preliminaries
% clc; clear all; close all;
%% Load in the data

function datamean = PLSDAFunction_June2022(filename)

dataPath = filename;

% addpath('C:\Working Copies\MATLAB Scripts\othercolor'); %Levi only.

% AMP 6/6/22 added to color correct images and error ellipse
addpath ('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\othercolor');
addpath ('error_ellipse');

% dataPath='20192002_Flicker240-263_cytokines_KMG_20190220_144041_inhibitors.xlsx';

[num,txt,raw]=xlsread(dataPath, 'Sheet1', 'D1:AI1');
analyteLabels=txt;

[num,txt,raw]=xlsread(dataPath, 'Sheet1', 'D4:AI27');
data=cell2mat(raw);

%Background wells
[num,txt,raw]=xlsread(dataPath, 'Sheet1', 'D2:AI3');
background=num;

[num,txt,raw]=xlsread(dataPath, 'Sheet1', 'B4:B27');
sampleLabels=raw;

[num,txt,raw]=xlsread(dataPath, 'Sheet1', 'C4:C27');
treatment=num;

%[num,txt,raw]=xlsread(dataPath, 'Sheet1', 'C2:C9');
%brainLoc=num;

%Subtract background
background=mean(background,1);
for i=1:length(background)  
    data(:,i)=data(:,i)-background(i);
end


%choose CTX Samples
%ind=find(brainLoc== 1)

%data=data(ind,:);
%sampleLabels=sampleLabels(ind);
%treatment=treatment(ind);
%brainLoc=brainLoc(ind);


%Set values less than zero to 0.
data(data<0)=0;



%z-score the data
zData=zscore(data);

%take a  look at the lysate panel for each brain region
% figure1= figure('units','normalized','position',[.1 .1 .45 .6])

%Pad the data, so that pcolor will work right
dataFig=zeros(size(data)+1); dataFig(1:end-1,1:end-1)=zscore(data);
pcolor(dataFig)
set(gca, 'ytick', [1:size(data,1)]+0.5)
set(gca, 'YTickLabel', sampleLabels)
set(gca, 'xtick', [1:length(analyteLabels)]+0.5)
set(gca, 'XTickLabel', analyteLabels)
set(gca,'position',[0.1 0.25 0.8 0.65])

th=rotateticklabel(gca, 90);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
c=othercolor('RdBu10');
colormap(flipud(c));
colorbar;
caxis([-2.5 2.5])


%PLSR
%//////////////////////////////////////////////////////////////////////


%pheno=treatment;
region.pls = pls(zData,treatment,3,'da');

%region.pls.P(:,1)=[]; region.pls.T(:,1)=[];


%Store loadings before the rotation
P1Store=region.pls.P(:,1);
P2Store=region.pls.P(:,2);

%Rotate the first two axes to see the seperating LV more clearly
%theta=45+90+25; % in degrees
theta=180; % in degrees

T1Temp=region.pls.T(:,1)*cos(pi/180*theta)+...
    region.pls.T(:,2)*sin(pi/180*theta);
T2Temp=region.pls.T(:,2)*cos(pi/180*theta)+...
    -region.pls.T(:,1)*sin(pi/180*theta);

region.pls.T(:,1)=T1Temp;
region.pls.T(:,2)=T2Temp;

%Rotate the loadings to match the scores
P1Temp=region.pls.P(:,1)*cos(pi/180*theta)+...
    region.pls.P(:,2)*sin(pi/180*theta);

P2Temp=region.pls.P(:,2)*cos(pi/180*theta)+...
    -region.pls.P(:,1)*sin(pi/180*theta);

region.pls.P(:,1)=P1Temp;
region.pls.P(:,2)=P2Temp;


ind40HzMAPKInh=find(treatment==1);
ind40HzVEH=find(treatment==2);
ind40HzNFkBInh=find(treatment==3);
ind20HzVEH=find(treatment==4);
%//////////////////////////////////////////////////////////////////////


%plot scores LV1-LV2
%//////////////////////////////////////////////////////////////////////
% fig2=figure('units','normalized','position',[.1 .2 [.2 .25]]);
% axes1 = axes('Parent',fig2,'LineWidth',2,'XTick',[-0.5 -0.25 0 0.25 0.5],'FontSize',13,'CLim',[0 3]);
% 
% plot(region.pls.T(ind40HzMAPKInh,1),region.pls.T(ind40HzMAPKInh,2),'o', 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]); hold on
% plot(region.pls.T(ind40HzVEH,1),region.pls.T(ind40HzVEH,2), 'o','MarkerFaceColor', [0 0 1], 'MarkerEdgeColor',[0 0 1])
% plot(region.pls.T(ind40HzNFkBInh,1),region.pls.T(ind40HzNFkBInh,2), 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
% plot(region.pls.T(ind20HzVEH,1),region.pls.T(ind20HzVEH,2), 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])


% xlim([-0.75 0.75])
% ylim([-0.75 0.75])
% box on
% xlabel(['Scores on LV1'], 'FontSize', 13);
% ylabel(['Scores on LV2'],'FontSize', 13)
% legend('40Hz MAPKInh','40Hz Veh','40Hz NFkBInh','20Hz Veh')
% set(gca, 'LineWidth', 1.5,'XTick',[-0.6:0.2:0.6])
% title('LV1 vs LV2')

%plot scores LV2-LV3
% fig3=figure('units','normalized','position',[.1 .2 [.2 .25]]);
% axes1 = axes('Parent',fig3,'LineWidth',2,'XTick',[-0.5 -0.25 0 0.25 0.5],'FontSize',13,'CLim',[0 3]);
% 
% plot(region.pls.T(ind40HzMAPKInh,2),region.pls.T(ind40HzMAPKInh,3),'o', 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]); hold on
% plot(region.pls.T(ind40HzVEH,2),region.pls.T(ind40HzVEH,3), 'o','MarkerFaceColor', [0 0 1], 'MarkerEdgeColor',[0 0 1])
% plot(region.pls.T(ind40HzNFkBInh,2),region.pls.T(ind40HzNFkBInh,3), 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
% plot(region.pls.T(ind20HzVEH,2),region.pls.T(ind20HzVEH,3), 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])

% xlim([-0.75 0.75])
% ylim([-0.75 0.75])
% box on
% xlabel(['Scores on LV2'], 'FontSize', 13);
% ylabel(['Scores on LV3'],'FontSize', 13)
% legend('40Hz MAPKInh','40Hz Veh','40Hz NFkBInh','20Hz Veh')
% set(gca, 'LineWidth', 1.5,'XTick',[-0.6:0.2:0.6])
% title('LV2 vs LV3')

%plot scores LV1-LV2-LV3
% fig4=figure('units','normalized','position',[.1 .2 [.2 .25]]);
% axes1 = axes('Parent',fig4,'LineWidth',2,'XTick',[-0.5 -0.25 0 0.25 0.5],'FontSize',13,'CLim',[0 3]);
% plot3(region.pls.T(ind40HzMAPKInh,1), region.pls.T(ind40HzMAPKInh,2),region.pls.T(ind40HzMAPKInh,3),'o', 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]); hold on
% plot3(region.pls.T(ind40HzVEH,1),region.pls.T(ind40HzVEH,2),region.pls.T(ind40HzVEH,3), 'o','MarkerFaceColor', [0 0 1], 'MarkerEdgeColor',[0 0 1])
% plot3(region.pls.T(ind40HzNFkBInh,1),region.pls.T(ind40HzNFkBInh,2),region.pls.T(ind40HzNFkBInh,3), 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
% plot3(region.pls.T(ind20HzVEH,1),region.pls.T(ind20HzVEH,2),region.pls.T(ind20HzVEH,3), 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])

% xlim([-0.75 0.75])
% ylim([-0.75 0.75])
% box on
% xlabel(['Scores on LV1'], 'FontSize', 13);
% ylabel(['Scores on LV2'],'FontSize', 13)
% zlabel(['Scores on LV3'],'FontSize', 13)
% legend('40Hz MAPKInh','40Hz Veh','40Hz NFkBInh','20Hz Veh')
% set(gca, 'LineWidth', 1.5,'XTick',[-0.6:0.2:0.6])
% title('LV1 vs LV2 vs LV3')


% bar plot the loadings on LV1
% fig5=figure('units','normalized','position',[.1 .2 [.35 .25]]);

% [region.pls.P(:,1) IX1]= sort(region.pls.P(:,1), 'ascend');
% bar([1:length(region.pls.P(:,1))], [region.pls.P(:,1)]/max(abs(region.pls.P(:,1))));%, ...
% set(gca, 'xtick', [1:length(analyteLabels)])
% set(gca, 'xTickLabel', analyteLabels (IX1))
% set(gca,'position',[0.1 0.35 0.8 0.55])
% h=gca;
% th=rotateticklabelLV(h, 90);
% title('Signals in LV1')
% xlim([0.5 size(data,2)+0.5]); ylim([-1 1])



% bar plot the loadings on LV2
% fig6=figure('units','normalized','position',[.1 .2 [.35 .25]]);

% [region.pls.P(:,2) IX2]= sort(region.pls.P(:,2), 'ascend');
% bar([1:length(region.pls.P(:,2))], [region.pls.P(:,2)]/max(abs(region.pls.P(:,2))));%, ...
% set(gca, 'xtick', [1:length(analyteLabels)])
% set(gca, 'xTickLabel', analyteLabels (IX2))
% set(gca,'position',[0.1 0.35 0.8 0.55])
% h=gca;
% th=rotateticklabelLV(h, 90);
% title('Signals in LV2')
% xlim([0.5 size(data,2)+0.5]); ylim([-1 1])


% clustergram(zData, 'RowLabels', sampleLabels, 'ColumnLabels', analyteLabels, 'ColorMap', flipud(c), 'Cluster', 'Row')

% figure('units','normalized','position',[.1 .1 .45 .6]);

%Pad the data, so that pcolor will work right
% dataFig=zeros(size(data)+1); dataFig(1:end-1,1:end-1)=zscore(data(:,IX1));
% pcolor(dataFig)
% set(gca, 'ytick', [1:size(data,1)]+0.5)
% set(gca, 'YTickLabel', sampleLabels)
% set(gca, 'xtick', [1:length(analyteLabels)]+0.5)
% set(gca, 'XTickLabel', analyteLabels(IX1))
% set(gca,'position',[0.1 0.25 0.8 0.65])

% th=rotateticklabel(gca, 90);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [8 8]);
% c=othercolor('RdBu10');
% colormap(flipud(c));
% colorbar;
% caxis([-2 2])

%% Replot cytokines in terms of LV1
% figure('units','normalized','position',[.1 .1 .45 .6]);

%Pad the data, so that pcolor will work right
% dataFig=zeros(size(data)+1); dataFig(1:end-1,1:end-1)=zscore(data(:,IX1));
% pcolor(dataFig)
% set(gca, 'ytick', [1:size(data,1)]+0.5)
% set(gca, 'YTickLabel', sampleLabels)
% set(gca, 'xtick', [1:length(analyteLabels)]+0.5)
% set(gca, 'XTickLabel', analyteLabels(IX1))
% set(gca,'position',[0.1 0.25 0.8 0.65])
% 
% th=rotateticklabel(gca, 90);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [8 8]);
% c=othercolor('RdBu10');
% colormap(flipud(c));
% colorbar;
% caxis([-2 2])


%use this for permutation averages per group AMP 6/27/22
 % separate LV scores for groups
 
%  [BarLabels,data,treatment]=unique(Group);
BarLabels = unique(treatment);
BarLabels=string(BarLabels);

% %  Table.Variable(Table.Variable == -1) = 0.1;
%  
 GroupLVscores={};
 datamean=[];
 for i=1:length(BarLabels);
     ind= find(treatment==i);
     GroupLVscores{i}=region.pls.T(ind,(1:2));
 end
%  
%  % get mean for each Group LV {'MAPKinhibitor+40Hz', 'VEH+40Hz', 'NFkBInhibitor+40Hz','VEH+20Hz'};

 datamean=[];
 for i=1:length(GroupLVscores);
     datamean(i,:)=mean(GroupLVscores{1,i});
 end


end

