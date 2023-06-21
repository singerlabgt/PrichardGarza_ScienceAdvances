%% Unblind data below run from Inhibitor Code Folder
% need this for process averages as opposed to direct data that doesn't
% need averages

% mousekey = readtable('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\key.xlsx');
mousekey = readtable('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\key.xlsx');
% mousekey = readtable('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\key.xlsx');

% string(mousekey.Mouse);
kh = height(mousekey); % how many mice in key file
uh = height(averages_table); % how much data to parse
placeholder=string(nan(uh,1)); % dummy array of Nans to populate (1 for condition, 2 for conditiona and blind ID)
averages_table.Subject=string(averages_table.Subject);

x = 1; % start iteration at 1

while(x<=kh)
    keymouse = string(mousekey.Mouse(x)); % for this mouse in list
    idx = ismember(lower(averages_table.Subject),lower(keymouse)); % find index of mosekey sample name that matches data sample
    num_matches = sum(idx); % number of matching rows from data
    matching_data_row_num=find(idx); % row index numbers from data that have mouse
    condition = string(mousekey.Condition(x));
    % create a line here if the mouse has a different blinded name
      
    placeholder(matching_data_row_num) = condition;
    x = x+1;
end

averages_table.Condition=placeholder; % updated with condition data

%% check for outliers based on condition
temp=table;
c = 1;
allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
while (c<=length(allConditions))
    group=allConditions(c);
    select_idx = averages_table.Condition==group ; % logical indexing 
    selection = averages_table(select_idx,:);
    selection.(5) = filloutliers(selection.(5),NaN); %replace outliers from median with NaNs
    temp=[temp; selection];
    c=c+1;
end
updatedinfo=temp;
updated_unblindedinfo=temp;
% updated_unblindedinfo=updatedinfo;
updatedfile=([logpathname 'Data\' region '_averages_' stringToBeFound '_unblinded_appended_' date '.mat']);
save(updatedfile,'updated_unblindedinfo');
disp(['saved unblinded data ' updatedfile]);
save([logpathname 'Data\updated_unblindedinfo.mat'],'updatedinfo');

%% averages stats

figure;
% Anova
[p,t,stats] = anovan(averages_table.(5), {averages_table.Condition}, "Varnames", "Condition");
manovafile = ([logpathname '\Stats\Anova_averages_' stringToBeFound '_' date ]);
save(manovafile,'t');

% run posthoc comparisons between groups 
posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
mmulitifile = ([logpathname '\Stats\TukeyKramerComps_averages_' stringToBeFound '_' date]);
save (mmulitifile, 'posthoc_comps');
disp('Anova and post-hoc comparison stats saved');
        
%% Means and Stdev
% [G,conditions] = findgroups(updatedinfo.Condition);
% avgdata = splitapply(@mean,updatedinfo.(4),G);
% stdevdata = splitapply(@std,updatedinfo.(4),G);
% avgdata.(1) = stringToBeFound;

[means,sems] = grpstats(averages_table.(5),{averages_table.Condition, averages_table.Subject},["mean","sem"]);
avgdata=table(means, sems);
gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
% avgdata.Properties.RowNames = gnames;
save(([logpathname '\Stats\Averages_' stringToBeFound '_' date]), 'avgdata');
disp('averages and sems saved');

%% create a violin plot for each value
disp('making graphs now');
    
% set path to save to
figpathname = ([pwd '\Graphs\']);
addpath(figpathname);

% add violin plot folder to path
addpath('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Violinplot-Matlab-master');
% addpath('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Violinplot-Matlab-master');
% addpath('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Violinplot-Matlab-master');

% "% plot the Violin", somewhere around line 145, you can add 'Linestyle', 'none' to the "fill" function and that will remove the outline
%% group Violin Plots
% figure;
%     averages_table.(3)=double(averages_table.(3));
%options ViolinColor, EdgeColor, MedianColor, BoxColor..
figure('Position', [19 402 462 398]);
vplot = violinplot(averages_table.(5), averages_table.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
vplot(1).ViolinColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(2).ViolinColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
vplot(2).EdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
vplot(3).ViolinColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
vplot(3).EdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
vplot(4).ViolinColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
vplot(4).EdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue

xlabel('Flicker Condition');
ylabel(['Average' stringToBeFound]);
saveas(gcf,[figpathname stringToBeFound '_averages_violin_plot_' date '.png']);

%% individual Violin Plots
% clf;
% figure;
allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
conditionColors = [0.6350 0.0780 0.1840 ; 0.9290 0.6940 0.1250 ;  0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330]; 

figure('Position', [-1554 512 1127 281]);
hold on
[subj, subjI] = unique(averages_table.Subject);
subjCond = averages_table.Condition(subjI);
subjColorMat = NaN(length(subj), 3);
for c = 1:length(allConditions)
    isCond = strcmp(subjCond, allConditions{c});
    subjColorMat(isCond,:) = repmat(conditionColors(c,:), sum(isCond),1);
end

vplot = violinplot(averages_table.(5), averages_table.Subject, 'BoxColor', [ 0 0 0]);
for c = 1:size(subjColorMat,1)
    vplot(c).ViolinColor = [subjColorMat(c,:)];
    vplot(c).EdgeColor = [subjColorMat(c,:)];
end

xlabel('Flicker Condition');
ylabel(stringToBeFound);
saveas(gcf,[figpathname stringToBeFound '_averages_individual_violin_plots_' date '.png']);
% end


%% Branch Depth/ Process Volume Ratio
load('VC_Filament_Full_Branch_Depth__appended_unblinded11-Apr-2023.mat', 'updated_unblindedinfo');
branchdepth_table=updated_unblindedinfo;
Dendrite_Volume_averages_table.mat;
load(branchdepth_table);
branch_by_volume=averages_table;
branch_by_volume.MeanProcessVolume=averages_table.MeanProcessVolume;
branch_by_volume.NumberofSegments=averages_table.NumberofSegments;
branch_by_volume.Ratio=branch_by_volume.Filament_Full_Branch_Depth ./ branch_by_volume.MeanProcessVolume;
save([logpathname 'Data\branch_by_volumeRatio_appended_04132023.mat'],'branch_by_volume');

% statistics
figure;
% Anova
[p,t,stats] = anovan(branch_by_volume.Ratio, {branch_by_volume.Condition}, "Varnames", "Condition");
manovafile = ([logpathname '\Stats\Anova_BranchDepthbyVolume_Ratio_' date ]);
save(manovafile,'t');

% run posthoc comparisons between groups 
posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
mmulitifile = ([logpathname '\Stats\TukeyKramerComps_BranchDepthbyVolume_Ratio_' date]);
save (mmulitifile, 'posthoc_comps');
disp('Anova and post-hoc comparison stats saved');

% [means,sems] = grpstats(branch_by_volume.Ratio,{branch_by_volume.Condition, branch_by_volume.Subject},["mean","sem"]);
avgdata=grpstats(branch_by_volume,["Condition","Subject"],["mean","sem"],"DataVars",'Ratio');

% avgdata=table(means, sems);
% gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
% avgdata.Properties.RowNames = gnames;
save(([logpathname '\Stats\Averages_BranchDepthbyVolume_Ratio_' date]), 'avgdata');
disp('averages and sems saved');

%% create a violin plot for each Ratio value
disp('making graphs now');
    
% group Violin Plots
% figure;
%     averages_table.(3)=double(averages_table.(3));
%options ViolinColor, EdgeColor, MedianColor, BoxColor..
figure('Position', [19 402 462 398]);
vplot = violinplot(branch_by_volume.Ratio, branch_by_volume.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
vplot(1).ViolinColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(2).ViolinColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
vplot(2).EdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
vplot(3).ViolinColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
vplot(3).EdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
vplot(4).ViolinColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
vplot(4).EdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
xlabel('Flicker Condition');
ylabel('Branching Depth/ Average Process Volume');
saveas(gcf,[figpathname 'Averages_BranchDepthbyVolume_Ratio_' date '.png']);

%% Branch Depth/ Soma Volume Ratio
load('VC_Soma_Volume__appended_unblinded11-Apr-2023.mat'), 'updated_unblindedinfo');
soma_volume_table=updated_unblindedinfo;
branch_by_volume.SomaVolume=soma_volume_averages_table.Soma_Volume;
branch_by_volume.SomaRatio=branch_by_volume.Filament_Full_Branch_Depth ./ branch_by_volume.SomaVolume;
save([logpathname 'Data\branch_by_SomavolumeRatio_appended_04132023.mat'],'branch_by_volume');

% statistics
figure;
% Anova
[p,t,stats] = anovan(branch_by_volume.SomaRatio, {branch_by_volume.Condition}, "Varnames", "Condition");
manovafile = ([logpathname '\Stats\Anova_BranchDepthby_SomaVolume_Ratio_' date ]);
save(manovafile,'t');

% run posthoc comparisons between groups 
posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
mmulitifile = ([logpathname '\Stats\TukeyKramerComps_BranchDepthby_SomaVolume_Ratio_' date]);
save (mmulitifile, 'posthoc_comps');
disp('Anova and post-hoc comparison stats saved');

% [means,sems] = grpstats(branch_by_volume.Ratio,{branch_by_volume.Condition, branch_by_volume.Subject},["mean","sem"]);
avgdata=grpstats(branch_by_volume,["Condition","Subject"],["mean","sem"],"DataVars",'SomaRatio');

% avgdata=table(means, sems);
% gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
% avgdata.Properties.RowNames = gnames;
save(([logpathname '\Stats\Averages_BranchDepthbySomaVolume_Ratio_' date]), 'avgdata');
disp('averages and sems saved');

%% create a violin plot for each Ratio value
disp('making graphs now');
    
% group Violin Plots
% figure;
%     averages_table.(3)=double(averages_table.(3));
%options ViolinColor, EdgeColor, MedianColor, BoxColor..
figure('Position', [19 402 462 398]);
vplot = violinplot(branch_by_volume.SomaRatio, branch_by_volume.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
vplot(1).ViolinColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(2).ViolinColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
vplot(2).EdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
vplot(3).ViolinColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
vplot(3).EdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
vplot(4).ViolinColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
vplot(4).EdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue

xlabel('Flicker Condition');
ylabel('Branching Depth/ Soma Volume');
saveas(gcf,[figpathname 'Averages_BranchDepthby_SomaVolume_Ratio_' date '.png']);

