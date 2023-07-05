%% create a violin plot for each value
clf; %clears current figure
disp('making graphs now');
% 
% % set path to save to
% figpathname = ('T:\singer\Ashley\ProcessedData\Figures\');
% addpath(figurepathname);

%% set up paths
addpath(genpath(pwd)); %code location
% Read in data
fileinfo = 'Flicker Colocalization Data Unblinded_with Cell Volume.xlsx';
% opts = detectImportOptions('Flicker Colocalization Data Unblinded_with Cell Volume.xlsx' , 'Sheet', 'Standardized Parameters');
% opts.SelectedVariableNames = ["ImageName", "TotalNFkBVolume", "CellVolume", "Condition"];
% T = readtable("patients.xls", opts)
data = readtable('Flicker Colocalization Data Unblinded_with Cell Volume.xlsx', 'Sheet', 'Standardized Parameters', 'ReadVariableNames', true, 'Range', 'A1:P46');
% data.Properties.VariableNames = 'TotalNFkBVolume';

% add violin plot folder to path
addpath('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Violinplot-Matlab-master');

%% remove blank data and groups we do not want to use below
toDelete = isnan(data.TotalNFkBVolume);
data(toDelete,:) = [];
novolume = isnan(data.CellVolume);
data(novolume,:) = [];

data(ismember(data.Condition,'40Hz PLX'),:)=[]; % remove PLX animals
data(ismember(data.Condition,'Light PLX'),:)=[]; % remove PLX animals
data(ismember(data.Condition,'Rand'),:)=[]; % remove PLX animals
% data(ismember(data.Condition,'Light'),:)=[]; % remove Light only control animals
data=data(:,[3 4 5 6 7 13 15]); % 15 is NFkB volume divided by cell volume

% remove 700's animals
% data(ismember(data.Subject,'701'),:)=[]; % remove extra animals
% data(ismember(data.Subject,'704'),:)=[]; % remove extra animals

data.Condition = categorical(data.Condition);
% location = {'NEUN'};
% data(ismember(data.Subject,'281'),:)=[]; % remove IBA1 outlier animals
% data(ismember(data.Subject,'283'),:)=[]; % remove IBA1 outlier animals

volume_outliers = data.nfkbColocalizedByCellVolume > 1.0; % remove outliers
data(volume_outliers,:) = [];

% location = {'NEUN'};
location = {'IBA1', 'NEUN'};

% count = 1;
for i = 1:length(location)
    loc = (location{i});
    x = ismember(data.Stain,loc);
    groupdata = data(x,:);

%% graph NFkB colocalization in violin plot
%     clf;
% %     clr = [1 0 0; 0 0 1; 0.5 0.5 0.5];
%     vplot = violinplot(groupdata.RatioColocalized, groupdata.Condition);
%     vplot(1).ViolinPlot.FaceColor = [0 0 1];
%     vplot(2).ViolinPlot.FaceColor = [1 0 0];
%     vplot(3).ViolinPlot.FaceColor = [0.5 0.5 0.5];
%     
%     vplot(1).ViolinPlot.MarkerFaceColor = [0 0 1];
%     vplot(2).ViolinPlot.MarkerFaceColor = [1 0 0];
%     vplot(3).ViolinPlot.MarkerFaceColor = [0.5 0.5 0.5];
%     
%     ylim([0 1.0]);
%     xlabel('Flicker Condition');
%     ylabel(['Percent NFkB Colocalized with ' loc ]);
    
%     saveas(gcf,['vplot_colocalized_NFkB_in_' loc '_' date '.png']);
%% Bar graphs option    
    clf;
    groups = unique(groupdata.Condition);
%     tblstats = grpstats(groupdata.RatioColocalized, groupdata.Condition); %group stats
    tblstats = grpstats(groupdata.nfkbColocalizedByCellVolume, groupdata.Condition);
    
%% calculate standard error of the mean for each group 
    rowslight = (groupdata.Condition=='Light'); % 40Hz rows
    LightHz=groupdata(rowslight,:);%get data from 40Hz
%     sem_40hz = std( FourHz.RatioColocalized) / sqrt (length( FourHz.RatioColocalized)); % calculate 40Hz SEM
    sem_light = std( LightHz.nfkbColocalizedByCellVolume) / sqrt (length( LightHz.nfkbColocalizedByCellVolume)); % calculate 40Hz SEM
%     bar_sem = [ sem_20hz sem_40hz];
    
    %%% stderror = std( data ) / sqrt( length );
    rows20 = (groupdata.Condition=='20Hz'); %20Hz rows
    TwHz=groupdata(rows20,:); %get data from 20Hz
%     sem_20hz = std( TwHz.RatioColocalized ) / sqrt( length( TwHz.RatioColocalized )); % calculate 20Hz SEM
    sem_20hz = std( TwHz.nfkbColocalizedByCellVolume ) / sqrt( length( TwHz.nfkbColocalizedByCellVolume )); % calculate 20Hz SEM

    rows40 = (groupdata.Condition=='40Hz'); % 40Hz rows
    FourHz=groupdata(rows40,:);%get data from 40Hz
%     sem_40hz = std( FourHz.RatioColocalized) / sqrt (length( FourHz.RatioColocalized)); % calculate 40Hz SEM
    sem_40hz = std( FourHz.nfkbColocalizedByCellVolume) / sqrt (length( FourHz.nfkbColocalizedByCellVolume)); % calculate 40Hz SEM
    bar_sem = [  sem_20hz sem_40hz sem_light];
    
    barfig = bar(groups, tblstats);
    barfig.FaceColor = 'flat';
    barfig.CData(1,:) = [0 0.4470 0.7410]; %20Hz = blue
%     barfig.CData(2,:) = [0.6350 0.0780 0.1840]; %40Hz = red 
    barfig.CData(2,:) = [1 0 0]; %40Hz = red 
    barfig.CData(3,:) = [.7 .7 .7];

    ylim([0 1]);
    
    hold on % hold to add error bars
    er = errorbar(groups,tblstats,bar_sem); % add SEM error bars for each group
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    
    % add scatter for individuals 
    scat = scatter(groupdata.Condition, groupdata.nfkbColocalizedByCellVolume);
    scat.MarkerEdgeColor = [0 0 0];
%     scat.MarkerFaceColor = [0 0.4470 0.7410; 1 0 0]; % if want colors

    hold off
    
    xlabel('Flicker Condition');
    ylabel(['Percent NFkB Volumen Colocalized with ' loc ' Volume']);
    saveas(gcf,['bar_colocalized_NFkB_in_' loc '_' date '.png']);
    disp('all graphs made and saved under Figures folder');
    
%% run statistics for group comparisons
% Anova
    disp('starting Anova and posthoc comps');

%     [p,t,stats] = anovan(groupdata.RatioColocalized, {groupdata.Condition}, "Varnames", "Condition");
    [p,t,stats] = anovan(groupdata.nfkbColocalizedByCellVolume, {groupdata.Condition}, "Varnames", "Condition");
    manovafile = (['Anova_NFkBcolocalization_' loc '_' date ]);
    save(manovafile,'t');
% 
% % run posthoc comparisons between groups for tukey HSD (standard setup)
    posthoccomps = multcompare(stats); %returns group, control group, lower limit, diff, upper limit, pval
    posthocfile = (['TukeyHSD_colocalization_' loc '_' date]);
    save (posthocfile, 'posthoccomps');
    disp('Anova and posthoc stats saved');

    % Fisher's exact test
%     G = groupsummary(data,{'Condition','Stain'},'mean','ColocalizedNFkBVolume');
%     z=table2array(G(:,4));
%     z=round(z);
%     G.mean_ColocalizedNFkBVolume = z;
%     contgtable(1,:) = G(1,4);
%     contgtable(2,:) = G(3,4);
%     contgtable(1,2) = G(2,4);
%     contgtable(2,2) = G(4,4);
%     contgtable.Properties.VariableNames={'IBA1', 'NEUN'};
%     contgtable.Properties.RowNames={'20Hz', '40Hz'};   
%     [h,p,stats] = fishertest(contgtable); % won't work for decimals
%     fishfile = (['FishersExact_NFkBcolocalization_' loc '_' date ]);
%     save(fishfile,'p');
    
    % Chi Squared cross tab
%     [conttbl,chi2,p,labels] = crosstab(data.Stain, data.Condition);
%     heatmap(data,'Stain','Condition');
%     chifile = (['Chi2_NFkBcolocalization_' date ]);
%     save(chifile,'p');
    
    i = i + 1; % increase the iteration
    disp(['done with ' loc]);
    disp(' ');
end

% attempt a linear mixed model - good to check if comparing regions
% must have subject IDs in dataset as they are repeated measures
% works (data, 'dependent variable ~ IV + (random effects | grouping variable)'
% lme = fitlme(data,'RatioColocalized ~ 1 + Stain + (1|Subject)');
% lme2 = fitlme(data,'RatioColocalized~Condition+(Condition|Stain)');
% 
% lme = fitlme(data,'RatioColocalized~Stain+(1|Condition)+(Stain-1|Condition)');
% lme = fitlme(data,'RatioColocalized~Stain+(1|Subject)+(Stain-1|Condition)');
% 
% lme2 = fitlme(data,'RatioColocalized~Stain+(Condition*Stain)');

disp('done');

% [p,t,stats] = anovan(data.RatioColocalized, {data.Stain}, "Varnames", "Stain");
% [p,t,stats] = anovan(data.RatioColocalized, {data.Condition}, "Varnames", "Condition");
%% overall anova results
% [p,t,stats] = anovan(data.RatioColocalized, {data.Stain, data.Condition}, "Varnames", {'Stain', 'Condition'});
% overall ANOVA across regions accounting for flicker conditon (20 or 40Hz)
% data(ismember(data.Condition,'Light'),:)=[]; % remove Light only control animals
[p,t,stats] = anovan(data.nfkbColocalizedByCellVolume, {data.Stain, data.Condition}, "Varnames", {'Stain', 'Condition'});
manovafile_overall = (['Anova_NFkBcolocalization_bothIBA1andNeuron_' date ]);
save(manovafile_overall,'t');

posthoccomps_overall = multcompare(stats); %returns group, control group, lower limit, diff, upper limit, pval
posthocfile = (['TukeyHSD_colocalization_overallcomps_' date]);
save (posthocfile, 'posthoccomps_overall');
disp('Anova and posthoc stats saved');


% tblstats = grpstats(data,"RatioColocalized","mean","DataVars",["Stain","Condition"]);
% tblstats = grpstats(data,"nfkbColocalizedByCellVolume","mean","DataVars",["Stain","Condition"]);
data.Stain = categorical(data.Stain);
GroupSum_conditionxstain = groupsummary(data,{'Condition','Stain'},'mean','nfkbColocalizedByCellVolume');
GroupSum_stainonly = groupsummary(data,{'Stain'},'mean','nfkbColocalizedByCellVolume');
rowsIBA1 = (data.Stain=='IBA1'); % IBA1 rows
rowsNEUN = (data.Stain=='NEUN'); % NeuN rows
dataIBA1=data(rowsIBA1,:);%get data from iba1
dataNEUN=data(rowsNEUN,:);%get data from neun
sem_IBA1 = std( dataIBA1.nfkbColocalizedByCellVolume) / sqrt (length( dataIBA1.nfkbColocalizedByCellVolume)); % calculate 40Hz SEM
sem_NEUN = std( dataNEUN.nfkbColocalizedByCellVolume) / sqrt (length( dataNEUN.nfkbColocalizedByCellVolume)); % calculate 40Hz SEM
