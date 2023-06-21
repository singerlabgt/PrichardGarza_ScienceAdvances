function []=graphing_function(updated_unblindedinfo,stringToBeFound,logpathname,updatedfile)
clf; %clears current figure
%% get paths run from Inhibitor Code Folder
addpath(genpath(pwd)); %code location
datapath = ([pwd '\Data\']);
addpath(datapath);
info=load('updated_unblindedinfo.mat');
updated_unblindedinfo=info.updatedinfo;

% set path to save to
figpathname = ([pwd '\Graphs\']);
addpath(figpathname);

%% Statistics
% clf; %clears current figure
figure;
% Anova
[p,t,stats] = anovan(updated_unblindedinfo.(4), {updated_unblindedinfo.Condition}, "Varnames", "Condition");
manovafile = ([logpathname '\Stats\Anova_' stringToBeFound '_' date ]);
save(manovafile,'t');

% run posthoc comparisons between groups 
posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
mmulitifile = ([logpathname '\Stats\TukeyKramerComps_' stringToBeFound '_' date]);
save (mmulitifile, 'posthoc_comps');
disp('Anova and post-hoc comparison stats saved');

%% Means and Stdev
% [G,conditions] = findgroups(updated_unblindedinfo.Condition);
% avgdata = splitapply(@mean,updatedinfo.(4),G);
% stdevdata = splitapply(@std,updatedinfo.(4),G);
% avgdata.(1) = stringToBeFound;

[means,sems] = grpstats(updated_unblindedinfo.(4),{updated_unblindedinfo.Condition},["mean","sem"]);
avgdata=table(means, sems);
gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
avgdata.Properties.RowNames = gnames;
save(([logpathname '\Stats\Averages' stringToBeFound '_' date]), 'avgdata');
disp('averages and sems saved');


% %% graph averages for subject and condition
% avgdata=grpstats(updated_unblindedinfo,["Condition","Subject"],["mean","sem"],"DataVars",stringToBeFound);
% % creage averages graph where each dot represents one animal
% figure('Position', [19 402 462 398]);
% vplot = violinplot(avgdata.(4), avgdata.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
% vplot(1).ViolinColor = {[0.6350 0.0780 0.1840]}; %vehicle + 40Hz = red
% %vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
% vplot(2).ViolinColor = {[0.9290 0.6940 0.1250]}; %vehicle + light = yellow
% %vplot(2).EdgeColor = {[0.9290 0.6940 0.1250]}; %vehicle + light = yellow
% vplot(3).ViolinColor = {[0.4660 0.6740 0.1880]}; %NFkB inhibitor + 40Hz = green
% %vplot(3).EdgeColor = {[0.4660 0.6740 0.1880]}; %NFkB inhibitor + 40Hz = green
% vplot(4).ViolinColor = {[0.3010 0.7450 0.9330]}; % NFkB inhibitor + light lbue
% %vplot(4).EdgeColor = {[0.3010 0.7450 0.9330]}; % NFkB inhibitor + light lbue
% 
% xlabel('Flicker Condition');
% ylabel(['Average' stringToBeFound]);
% saveas(gcf,[figpathname stringToBeFound '_averages_violin_plot_' date '.png']);
% 
% figure;
% % Anova
% [p,t,stats] = anovan(avgdata.(4), {avgdata.Condition}, "Varnames", "Condition");
% manovafile = ([logpathname '\Stats\Anova_averages_' stringToBeFound '_' date ]);
% save(manovafile,'t');
% 
% % run posthoc comparisons between groups 
% posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
% mmulitifile = ([logpathname '\Stats\TukeyKramerComps_averages_' stringToBeFound '_' date]);
% save (mmulitifile, 'posthoc_comps');
% disp('Anova and post-hoc comparison stats saved');

%% graph averages for subject and condition
avgdata=grpstats(updated_unblindedinfo,["Condition","Subject","FilamentID"],["mean","sem"],"DataVars",stringToBeFound);
% creage averages graph where each dot represents one animal
figure('Position', [19 402 462 398]);
vplot = violinplot(avgdata.(4), avgdata.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
vplot(1).ViolinColor = {[0.6350 0.0780 0.1840]}; %vehicle + 40Hz = red
%vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(2).ViolinColor = {[0.9290 0.6940 0.1250]}; %vehicle + light = yellow
%vplot(2).EdgeColor = {[0.9290 0.6940 0.1250]}; %vehicle + light = yellow
vplot(3).ViolinColor = {[0.4660 0.6740 0.1880]}; %NFkB inhibitor + 40Hz = green
%vplot(3).EdgeColor = {[0.4660 0.6740 0.1880]}; %NFkB inhibitor + 40Hz = green
vplot(4).ViolinColor = {[0.3010 0.7450 0.9330]}; % NFkB inhibitor + light lbue
%vplot(4).EdgeColor = {[0.3010 0.7450 0.9330]}; % NFkB inhibitor + light lbue

xlabel('Flicker Condition');
ylabel(['Average' stringToBeFound]);
saveas(gcf,[figpathname stringToBeFound '_averages_violin_plot_' date '.png']);

figure;
% Anova
[p,t,stats] = anovan(avgdata.(4), {avgdata.Condition}, "Varnames", "Condition");
manovafile = ([logpathname '\Stats\Anova_averages_' stringToBeFound '_' date ]);
save(manovafile,'t');

% run posthoc comparisons between groups 
posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
mmulitifile = ([logpathname '\Stats\TukeyKramerComps_averages_' stringToBeFound '_' date]);
save (mmulitifile, 'posthoc_comps');
disp('Anova and post-hoc comparison stats saved');

%% create a violin plot for each value
disp('making graphs now');



% "% plot the Violin", somewhere around line 145, you can add 'Linestyle', 'none' to the "fill" function and that will remove the outline
%% group Violin Plots
% figure;
updated_unblindedinfo.(3)=double(updated_unblindedinfo.(3));
%options ViolinColor, EdgeColor, MedianColor, BoxColor..
figure('Position', [19 402 462 398]);
vplot = violinplot(updated_unblindedinfo.(4), updated_unblindedinfo.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
vplot(1).ViolinColor = {[0.6350 0.0780 0.1840]}; %vehicle + 40Hz = red
% vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
vplot(2).ViolinColor ={[0.9290 0.6940 0.1250]}; %vehicle + light = yellow
% vplot(2).EdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
vplot(3).ViolinColor = {[0.4660 0.6740 0.1880]}; %NFkB inhibitor + 40Hz = green
% vplot(3).EdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
vplot(4).ViolinColor = {[0.3010 0.7450 0.9330]}; % NFkB inhibitor + light lbue
% vplot(4).EdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue

xlabel('Flicker Condition');
ylabel(stringToBeFound);
saveas(gcf,[figpathname stringToBeFound '_violin_plot_' date '.png']);

%% individual Violin Plots
% clf;
% figure;
% allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
% conditionColors = {[0.6350 0.0780 0.1840 ; 0.9290 0.6940 0.1250 ;  0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330]}; 
% 
% figure('Position', [-1554 512 1127 281]);
% hold on
% [subj, subjI] = unique(updated_unblindedinfo.Subject);
% subjCond = updated_unblindedinfo.Condition(subjI);
% subjColorMat = NaN(length(subj), 3);
% for c = 1:length(allConditions)
%     isCond = strcmp(subjCond, allConditions{c});
%     subjColorMat(isCond,:) = repmat(conditionColors(c,:), sum(isCond),1);
% end
% 
% vplot = violinplot(updated_unblindedinfo.(4), updated_unblindedinfo.Subject, 'BoxColor', [ 0 0 0]);
% for c = 1:size(subjColorMat,1)
%     vplot(c).ViolinColor = [subjColorMat(c,:)];
%     vplot(c).EdgeColor = [subjColorMat(c,:)];
% end
% 
% xlabel('Flicker Condition');
% % stringToBeFound='Surface Area';
% ylabel(stringToBeFound);
% saveas(gcf,[figpathname stringToBeFound '_individual_violin_plots_' date '.png']);






%% Statistics
% clf; %clears current figure
% figure;
% % Anova
% [p,t,stats] = anovan(updated_unblindedinfo.(4), {updated_unblindedinfo.Condition}, "Varnames", "Condition");
% manovafile = ([logpathname '\Stats\Anova_' stringToBeFound '_' date ]);
% save(manovafile,'t');
% 
% % run posthoc comparisons between groups 
% posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
% mmulitifile = ([logpathname '\Stats\TukeyKramerComps_' stringToBeFound '_' date]);
% save (mmulitifile, 'posthoc_comps');
% disp('Anova and post-hoc comparison stats saved');

%% Means and Stdev
% [G,conditions] = findgroups(updatedinfo.Condition);
% avgdata = splitapply(@mean,updatedinfo.(4),G);
% stdevdata = splitapply(@std,updatedinfo.(4),G);
% avgdata.(1) = stringToBeFound;
% 
% [means,sems] = grpstats(updated_unblindedinfo.(4),updated_unblindedinfo.Condition,["mean","sem"]);
% avgdata=table(means, sems);
% gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
% avgdata.Properties.RowNames = gnames;
% save(([logpathname '\Stats\Averages' stringToBeFound '_' date]), 'avgdata');
% disp('averages and sems saved');
hold on
end