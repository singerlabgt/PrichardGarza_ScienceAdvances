%GetImarisData
%Code to Retrieve and Colosolidate Morphology Data from Imaris CV files
clc;
clear;
close all;
%% set up paths
addpath(genpath(pwd)); %code location run from Inhibitor Code Folder

% directory location for final data in Shared Box
% logpathname = ('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
logpathname = ('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
% logpathname = ('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');

% location of data in Emma's Neurocloud Folder
% datapathname = ('T:\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
datapathname = ('Y:\singer\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
addpath(genpath(datapathname)); %data location including all folders and subfolders for each animal

% add violin plot folder to path
violinplot_location=('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_code\Violinplot-Matlab-master');
addpath(violinplot_location);

% location of path to save new data to
newdatapath = ([pwd '\Data\']);
addpath(newdatapath);

% set path to save figures to
figpathname = ([pwd '\Graphs\']);
addpath(figpathname);

%% set up the data matrix
% create a new info table to populate with morphology values
originfo = array2table(NaN(1,4));
originfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'VariableofInterest'};
% use this line below if you have multiple samples for each mouse
% originfo = array2table(NaN(1,7));
% originfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};

%% edit below for region
region = 'VC';

%% morphology specifics to look for in the Imaris file name

% Area Measures (total cell area, soma area, average branch length) 
% microglia_number = 'Soma Area ImageJ.xlsx';
% soma_area = 'Soma Area Image J.xlsx'; %this was done in ImageJ don't use
% Imaris files

total_process_area = 'Filament_Area_(sum)'; % does not include soma
branching_depth = 'Filament_Full_Branch_Depth';
total_process_length = 'Filament_Length_(sum)'; % length per process
% for total_number_processes we can get number per cell

% Volume Measures (total cell volume, soma volume, average branch volume)
total_process_volume = 'Filament_Volume_(sum)'; % total process volume per microglia - the soma (no soma here)
soma_volume = 'Soma_Volume';
branch_points = 'Filament_No._Dendrite_Branch_Pts'; % number of nodes
number_of_branches = 'Filament_No._Dendrite_Branches'; %number of branches resulting from nodes

% segment measures (node to node segments and not the entire process)
average_segment_length = 'Dendrite_Length';
average_segment_volume = 'Dendrite_Volume';

% convex hull area
convexhull_area = 'hull_Area';

% convex hull volume
convexhull_volume = 'hull_Volume';


%% **** specify string to find here **** run one at a time
stringToBeFound = branch_points;
disp(['now calculating ' stringToBeFound ' for each sample in the ' region ' using means and maxes for each Filament ID...']);

%% read in folders and begin iteration through folder
filesAndfolders = dir(datapathname); % returns all files and folders in directory
numoffolders = length(filesAndfolders); % the number of folders in the directory
b = 1:size(filesAndfolders);
i=3; % set a folder to start, first 2 folders are filler '..' and should be skipped for this process

averages_table=table; % create table to update each iteration of averages
updatedinfo = originfo; %create table to update each iteration of raw values

while(i<=numoffolders)
    curdir = filesAndfolders(i).name;    % Store the name of the folder
%     mouseID = extractBetween(curdir,"(",")"); % get animal ID
%     mouseID=string(mouseID);
    mouseID=string(curdir);
    disp(['finding morphology data for' mouseID]);
%     if contains(region,'VC') % -- only need this if multiple samples per animal
%         sampleID= extractBetween(curdir,"40x_","CA1"); %get specimen number
%         sampleID=string(sampleID);
%         sampleID= extractBetween(sampleID,"_","_"); 
%     end
    mouseinfo = ([datapathname curdir]); %create a path to the folder
    mousedir = dir(mouseinfo); % make it a directory
    num_dir=numel(mousedir); % count number of files
    all_filesandfoldersindir = dir(fullfile(mouseinfo,'\**\*.*'));
    filesindir = all_filesandfoldersindir(~([all_filesandfoldersindir.isdir])); % not folders
    numberffiles = length(filesindir); % the number of files in the directory

    a=1; %create an interation
    while(a<=numberffiles) %start the while loop through files within the directory
        filename = filesindir(a).name;  % Store the name of the file
        found = strfind(filename, stringToBeFound); % look for specific variable in file name listed above
        if ~isempty(found) % if the results are not empty, it has found the file name
            foundString = strcat('Found in file ------', filename);
            disp(foundString); % sanity check

            % find data from first column
%             thismouse = ([mouseinfo '\' filename]);
            thismouse = filename;
            % read in excel key with subject data
            data = readtable(thismouse);
%             data.(1) = filloutliers(data.(1),NaN); %replace outliers from median with NaNs
            num_cells = height(data); % number of filament rows = cell count
            dummy_mouse = repmat(mouseID,num_cells,1); % populate mouse name
            dummy_region = repmat(string(region),num_cells,1); % populate cell data
            newinfo = table(dummy_mouse, dummy_region, data.ID, data.(1));
            newinfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'VariableofInterest' }; % filament ID = cell ID
            updatedinfo = [updatedinfo; newinfo]; % append rows to table for each subset of mouse data

            %% use below to get means by cell for each subject
            avgdata=grpstats(newinfo,["Subject","FilamentID"],["mean","sem"],"DataVars",'VariableofInterest');
            averages_table = [averages_table; avgdata];

            break;
        end
        a=a+1; % iterate through additional files in folder if not found
    end
    i = i+1; % move onto the next mouse ID
    
end        
%%        
disp(['updating info to table for ' stringToBeFound]); % display sanity check for mouse
% disp(['in the region ' region]); %display sanity check for region
disp('done');

updatedinfo(1,:) = []; %delete first blank row
avgcalc= strcmp(stringToBeFound,{'Dendrite_Length' , 'Dendrite_Volume'});
if sum(avgcalc) == 1
    avgdata_bydendrite=grpstats(updatedinfo,["Subject","FilamentID"],["mean","sem"],"DataVars",'VariableofInterest');
end

updatedinfo=renamevars(updatedinfo,"VariableofInterest",stringToBeFound); %tename variable to string of interest

updatedrawfile=([logpathname 'Data\' region '_' stringToBeFound '__appended_' date '.mat']); 
save(updatedrawfile,'updatedinfo'); % save raw values file

averages_table=renamevars(averages_table,"mean_VariableofInterest",['mean_' stringToBeFound]); % rename averages variable to string of interest

updatedaveragefile = ([logpathname 'Data\' region '_' stringToBeFound '_averages_appended_' date '.mat']); 
save(updatedaveragefile, 'averages_table'); % save average values file

disp('all data saved ');

%% Unblinding Code Below
disp('starting unblinding code below');  
% read in subject data from key
mousekey = readtable('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\key.xlsx');
number_of_mice = height(mousekey); % how many mice in key file
updated_data_height = height(updatedinfo); % how much data to parse

placeholder=string(nan(updated_data_height,1)); % dummy array of Nans to populate (1 for condition, 2 for conditiona and blind ID)
updatedinfo.Subject=string(updatedinfo.Subject); % make sure subject names are strings
%% create an iteration to add condition column to raw values table
x = 1; % start iteration at 1
while(x<=number_of_mice)
    keymouse = string(mousekey.Mouse(x)); % for this mouse in list
    idx = ismember(lower(updatedinfo.Subject),lower(keymouse)); % find index of mosekey sample name that matches data sample
    num_matches = sum(idx); % number of matching rows from data
    matching_data_row_num=find(idx); % row index numbers from data that have mouse
    condition = string(mousekey.Condition(x));
    % create a line here if the mouse has a different blinded name
      
    placeholder(matching_data_row_num) = condition;
    x = x+1;
end
updatedinfo.Condition=placeholder; % updated table column with condition data

%% create an iteration to add condition column to average values table
keyheight = height(mousekey); % how many mice in key file
averagedata_height = height(averages_table); % how much data to parse
placeholder=string(nan(averagedata_height,1)); % dummy array of Nans to populate (1 for condition, 2 for conditiona and blind ID)
averages_table.Subject=string(averages_table.Subject); %make sure subject names are strings

x = 1; % start iteration at 1

while(x<=keyheight)
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
disp('data unblinded');
%% Check for outliers based on condition in raw values table

rawtemp=table;
c = 1;
allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
while (c<=length(allConditions))
    group=allConditions(c);
    select_idx = updatedinfo.Condition==group ; % logical indexing 
    selection = updatedinfo(select_idx,:);
    selection.(4) = filloutliers(selection.(4),NaN); %replace outliers from median with NaNs
    rawtemp=[rawtemp; selection];
    c=c+1;
end
unblinded_rawvalues_table=rawtemp;
% updated_unblindedinfo=updatedinfo;
updatedfile=([logpathname 'Data\' region '_' stringToBeFound '__appended_unblinded' date '.mat']);
save(updatedfile,'unblinded_rawvalues_table');
disp('saved unblinded data');

%% check for outliers based on condition in averages value table
avgtemp=table;
c = 1;
allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
while (c<=length(allConditions))
    group=allConditions(c);
    select_idx = averages_table.Condition==group ; % logical indexing 
    selection = averages_table(select_idx,:);
    selection.(5) = filloutliers(selection.(5),NaN); %replace outliers from median with NaNs
    avgtemp=[avgtemp; selection];
    c=c+1;
end
unblinded_averages_table=avgtemp;
updated_avgfile=([logpathname 'Data\' region '_averages_' stringToBeFound '_unblinded_appended_' date '.mat']);
save(updated_avgfile,'unblinded_averages_table');
disp(['saved unblinded data ' updated_avgfile]);

%% Statistics
disp ('running statistics now')

% Anova for between groups comparison
[p,t,stats] = anovan(unblinded_rawvalues_table.(4), {unblinded_rawvalues_table.Condition}, "Varnames", "Condition");
manovafile = ([logpathname '\Stats\Anova_' stringToBeFound '_' date ]);
save(manovafile,'t');

% run Tukey's posthoc comparisons between groups 
posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
mmulitifile = ([logpathname '\Stats\TukeyKramerComps_' stringToBeFound '_' date]);
save (mmulitifile, 'posthoc_comps');
disp('Anova and post-hoc comparison stats saved');

%% Graphs
disp('making graphs now');
% "% plot the Violin", somewhere around line 145, you can add 'Linestyle', 'none' to the "fill" function and that will remove the outline
% create a violin plot for each condition based on data by cell

figure('Position', [19 402 462 398]);
vplot = violinplot(unblinded_rawvalues_table.(4), unblinded_rawvalues_table.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
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

%% create a violin plot for each individual within each condition by cell
% allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
% conditionColors = [0.6350 0.0780 0.1840 ; 0.9290 0.6940 0.1250 ;  0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330]; 
% 
% figure('Position', [-1554 512 1127 281]);
% hold on
% [subj, subjI] = unique(unblinded_rawvalues_table.Subject);
% subjCond = unblinded_rawvalues_table.Condition(subjI);
% subjColorMat = NaN(length(subj), 3);
% for c = 1:length(allConditions)
%     isCond = strcmp(subjCond, allConditions{c});
%     subjColorMat(isCond,:) = repmat(conditionColors(c,:), sum(isCond),1);
% end
% 
% vplot = violinplot(unblinded_rawvalues_table.(4), unblinded_rawvalues_table.Subject, 'BoxColor', [ 0 0 0]);
% for c = 1:size(subjColorMat,1)
%     vplot(c).ViolinColor = [subjColorMat(c,:)];
%     vplot(c).EdgeColor = [subjColorMat(c,:)];
% end
% 
% xlabel('Flicker Condition');
% ylabel(stringToBeFound);
% saveas(gcf,[figpathname stringToBeFound '_individual_violin_plots_' date '.png']);

%% if segment data is being used, also create a graph using averages
if sum(avgcalc) == 1
    % get averages
%     avgdata_bydendrite=grpstats(unblinded_rawvalues_table,["Subject","FilamentID"],["mean","sem"],"DataVars",'VariableofInterest');
    data_to_avg = averages_table.(4);
    avgdata_dendritebycell=grpstats(averages_table,["Condition","Subject","FilamentID"],["mean","sem"], "DataVars", 4);

    figure('Position', [19 402 462 398]);
    vplot = violinplot(averages_table.(4), averages_table.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
    vplot(1).ViolinColor = {[0.6350 0.0780 0.1840]}; %vehicle + 40Hz = red
%     vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
    vplot(2).ViolinColor ={[0.9290 0.6940 0.1250]}; %vehicle + light = yellow
%     vplot(2).EdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
    vplot(3).ViolinColor = {[0.4660 0.6740 0.1880]}; %NFkB inhibitor + 40Hz = green
%     vplot(3).EdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
    vplot(4).ViolinColor = {[0.3010 0.7450 0.9330]}; % NFkB inhibitor + light lbue
%     vplot(4).EdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
    ylim([0,(max(averages_table.(4))+1)]);
    xlabel('Flicker Condition');
    ylabel(['Average' stringToBeFound ' by Cell']);
    saveas(gcf,[figpathname stringToBeFound '_averages_violin_plot_' date '.png']);
end

% create segment histograms for Annabelle and Levi 4/15
if sum(avgcalc) == 1
    allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
%     conditionColors = [0.6350 0.0780 0.1840 ; 0.9290 0.6940 0.1250 ;  0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330]; 
    figure;
    c1= unblinded_rawvalues_table.Condition==allConditions(1);
    h1=histogram((unblinded_rawvalues_table.(4)(c1)));
    hold on
    c2= unblinded_rawvalues_table.Condition==allConditions(2);
    h2=histogram((unblinded_rawvalues_table.(4)(c2)));
    c3= unblinded_rawvalues_table.Condition==allConditions(3);
    h3=histogram((unblinded_rawvalues_table.(4)(c3)));
    c4= unblinded_rawvalues_table.Condition==allConditions(4);
    h4=histogram((unblinded_rawvalues_table.(4)(c4)));
    legend(allConditions);
    ylabel([stringToBeFound ' Count']);
    hold off;
%     h1.Normalization = 'probability';
%     h1.BinWidth = 0.25;
%     h2.Normalization = 'probability';
%     h2.BinWidth = 0.25;
%     h3.Normalization = 'probability';
%     h3.BinWidth = 0.25;
%     h4.Normalization = 'probability';
%     h4.BinWidth = 0.25;

 
    


%% get histogram in line form 4/16
    figure;
    [N1,edges] = histcounts((unblinded_rawvalues_table.(4)(c1)));
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N1);
    hold on
    [N2,edges] = histcounts((unblinded_rawvalues_table.(4)(c2)));
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N2);
    [N3,edges] = histcounts((unblinded_rawvalues_table.(4)(c3)));
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N3);
    [N4,edges] = histcounts((unblinded_rawvalues_table.(4)(c4)));
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N4);
    legend(allConditions)
    ylabel([stringToBeFound ' Count']);
    hold off
end   
    
%      for z = 1:length(allConditions)
%          Condition=allConditions(z);
%          thisdata = unblinded_rawvalues_table.Condition==Condition;
%          h1=histogram((unblinded_rawvalues_table.(4)(thisdata)));
%          hold on
%      end

%     isCond = strcmp(subjCond, allConditions{c});
%     subjColorMat(isCond,:) = repmat(conditionColors(c,:), sum(isCond),1);
% end

%% Branch Depth/ Process Volume Ratio
% bd=load('VC_Filament_Full_Branch_Depth__unblinded_appended_16-Apr-2023.mat', 'unblinded_rawvalues_table');
% branchdepth_table=unblinded_rawvalues_table;
% pv=load('VC_Filament_Volume_(sum)_unblinded_appended_16-Apr-2023.mat', 'unblinded_rawvalues_table');
% 
% % Dendrite_Volume_averages_table.mat;
% load(branchdepth_table);
% branch_by_volume=averages_table;
% branch_by_volume.MeanProcessVolume=averages_table.MeanProcessVolume;
% branch_by_volume.NumberofSegments=averages_table.NumberofSegments;
% branch_by_volume.Ratio=branch_by_volume.Filament_Full_Branch_Depth ./ branch_by_volume.MeanProcessVolume;
% save([logpathname 'Data\branch_by_volumeRatio_appended_04132023.mat'],'branch_by_volume');
% 
% % statistics
% figure;
% % Anova
% [p,t,stats] = anovan(branch_by_volume.Ratio, {branch_by_volume.Condition}, "Varnames", "Condition");
% manovafile = ([logpathname '\Stats\Anova_BranchDepthbyVolume_Ratio_' date ]);
% save(manovafile,'t');
% 
% % run posthoc comparisons between groups 
% posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
% mmulitifile = ([logpathname '\Stats\TukeyKramerComps_BranchDepthbyVolume_Ratio_' date]);
% save (mmulitifile, 'posthoc_comps');
% disp('Anova and post-hoc comparison stats saved');
% 
% % [means,sems] = grpstats(branch_by_volume.Ratio,{branch_by_volume.Condition, branch_by_volume.Subject},["mean","sem"]);
% avgdata=grpstats(branch_by_volume,["Condition","Subject"],["mean","sem"],"DataVars",'Ratio');
% 
% % avgdata=table(means, sems);
% % gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
% % avgdata.Properties.RowNames = gnames;
% save(([logpathname '\Stats\Averages_BranchDepthbyVolume_Ratio_' date]), 'avgdata');
% disp('averages and sems saved');
% 
% %% create a violin plot for each Ratio value
% disp('making graphs now');
%     
% % group Violin Plots
% % figure;
% %     averages_table.(3)=double(averages_table.(3));
% %options ViolinColor, EdgeColor, MedianColor, BoxColor..
% figure('Position', [19 402 462 398]);
% vplot = violinplot(branch_by_volume.Ratio, branch_by_volume.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
% vplot(1).ViolinColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
% vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
% vplot(2).ViolinColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
% vplot(2).EdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
% vplot(3).ViolinColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
% vplot(3).EdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
% vplot(4).ViolinColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
% vplot(4).EdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
% xlabel('Flicker Condition');
% ylabel('Branching Depth/ Average Process Volume');
% saveas(gcf,[figpathname 'Averages_BranchDepthbyVolume_Ratio_' date '.png']);
% 
% %% Branch Depth/ Soma Volume Ratio
% load('VC_Soma_Volume__appended_unblinded11-Apr-2023.mat'), 'updated_unblindedinfo');
% soma_volume_table=updated_unblindedinfo;
% branch_by_volume.SomaVolume=soma_volume_averages_table.Soma_Volume;
% branch_by_volume.SomaRatio=branch_by_volume.Filament_Full_Branch_Depth ./ branch_by_volume.SomaVolume;
% save([logpathname 'Data\branch_by_SomavolumeRatio_appended_04132023.mat'],'branch_by_volume');
% 
% % statistics
% figure;
% % Anova
% [p,t,stats] = anovan(branch_by_volume.SomaRatio, {branch_by_volume.Condition}, "Varnames", "Condition");
% manovafile = ([logpathname '\Stats\Anova_BranchDepthby_SomaVolume_Ratio_' date ]);
% save(manovafile,'t');
% 
% % run posthoc comparisons between groups 
% posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
% mmulitifile = ([logpathname '\Stats\TukeyKramerComps_BranchDepthby_SomaVolume_Ratio_' date]);
% save (mmulitifile, 'posthoc_comps');
% disp('Anova and post-hoc comparison stats saved');
% 
% % [means,sems] = grpstats(branch_by_volume.Ratio,{branch_by_volume.Condition, branch_by_volume.Subject},["mean","sem"]);
% avgdata=grpstats(branch_by_volume,["Condition","Subject"],["mean","sem"],"DataVars",'SomaRatio');
% 
% % avgdata=table(means, sems);
% % gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
% % avgdata.Properties.RowNames = gnames;
% save(([logpathname '\Stats\Averages_BranchDepthbySomaVolume_Ratio_' date]), 'avgdata');
% disp('averages and sems saved');
% 
% %% create a violin plot for each Ratio value
% disp('making graphs now');
%     
% % group Violin Plots
% % figure;
% %     averages_table.(3)=double(averages_table.(3));
% %options ViolinColor, EdgeColor, MedianColor, BoxColor..
% figure('Position', [19 402 462 398]);
% vplot = violinplot(branch_by_volume.SomaRatio, branch_by_volume.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light',  'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
% vplot(1).ViolinColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
% vplot(1).EdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
% vplot(2).ViolinColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
% vplot(2).EdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
% vplot(3).ViolinColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
% vplot(3).EdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
% vplot(4).ViolinColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
% vplot(4).EdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
% 
% xlabel('Flicker Condition');
% ylabel('Branching Depth/ Soma Volume');
% saveas(gcf,[figpathname 'Averages_BranchDepthby_SomaVolume_Ratio_' date '.png']);


