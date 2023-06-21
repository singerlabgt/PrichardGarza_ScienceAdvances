%Get ImageJ Data
%Code to Retrieve and Colosolidate Morphology Data from Imaris CV files
clc;
clear;

%% set up paths
 addpath(genpath(pwd)); %code location
%directory location for final data in Box
% logpathname = ('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
logpathname = ('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
% logpathname = ('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');

%location of data in Emma's Neurocloud Folder
datapathname = ('Y:\singer\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
% datapathname = ('T:\singer\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
addpath(genpath(datapathname)); %data location including all folders and subfolders for each animal

% create a new info table to populate
originfo = array2table(NaN(1,3));
originfo.Properties.VariableNames = {'Subject' 'Region' 'Soma Area'};
% use this line below if you have multiple samples for each mouse
% originfo = array2table(NaN(1,7));
% originfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};

%% edit below for region
region = 'VC';

%% morphology specifics to look for in the Imaris file name

soma_area = 'Soma Area';
% microglia_number = 'Soma Area ImageJ.xlsx';

%% **** specify string to find here **** run one at a time
stringToBeFound = soma_area;
disp(['now calculating ' stringToBeFound ' for each sample in the ' region ' using means and maxes for each Filament ID...']);

data=readtable([logpathname 'Data\Soma Area Image J.xlsx']);

%% read in column of data for each mouse
updatedinfo = originfo; %create table to update each iteration
number_of_mice = size(data,2);

number_of_glia=originfo;
number_of_glia=renamevars(number_of_glia,"Soma Area",'Number of Microglia');

i=1; % set a mouse to start, iterate by column

while(i<=number_of_mice)
    mouseID = data.Properties.VariableNames(i); %get the current column name
    mouseID=string(mouseID);
    curdir=data(:,i);%  % get the current column data
    curdir=rmmissing(curdir); % remove nans
    num_cells=height(curdir); % get number of cells (number of rows)
    curdir.Properties.VariableNames={'Soma Area'};

    placeholder_names=table(repmat(mouseID,num_cells,1)); % create matching number of names
    placeholder_names.Properties.VariableNames={'Subject'};
    placeholder_region=table(repmat(region,num_cells,1)); % create matching number for region
    placeholder_region.Properties.VariableNames={'Region'};
    tempinfo = [placeholder_names placeholder_region curdir]; % create temp table for new data
    tempinfo.Region=string(tempinfo.Region);       
    updatedinfo = [updatedinfo; tempinfo]; %append rows to table for each subset
    loc=string(region);        
    %% cell info
    celldata=table(mouseID, loc, num_cells);
    celldata.Properties.VariableNames={'Subject' 'Region' 'Number of Microglia'};
    number_of_glia=[number_of_glia; celldata];
    
    i = i+1; % move onto the next mouse ID
    
end        
%%        
disp(['updating info to table for ' stringToBeFound]);
disp(['in the region ' region]);
disp('done');

updatedinfo(1,:) = []; %delete first blank row
number_of_glia(1,:) = []; %delete first blank row

% updatedinfo.Filament_Full_Branch_Depth=updatedinfo.VariableofInterest;

updatedfile=([logpathname 'Data\' region '_' stringToBeFound '__appended_' date '.mat']);
save(updatedfile,'updatedinfo');

updatedfile=([logpathname 'Data\' region '_ number_of_glia__appended_' date '.mat']);
save(updatedfile,'number_of_glia');

save([logpathname 'Data\updatedinfo.mat'],'updatedinfo');
disp(['saved ' updatedfile]);
disp('passing to unblinding function');  

%% unblinding_function(updatedinfo,region,stringToBeFound,logpathname)
% mousekey = readtable('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\key.xlsx');
mousekey = readtable('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\key.xlsx');

% string(mousekey.Mouse);
kh = height(mousekey); % how many mice in key file
uh = height(updatedinfo); % how much data to parse
placeholder=string(nan(uh,1)); % dummy array of Nans to populate (1 for condition, 2 for conditiona and blind ID)
updatedinfo.Subject=string(updatedinfo.Subject);

x = 1; % start iteration at 1

while(x<=kh)
    keymouse = string(mousekey.Mouse(x)); % for this mouse in list
    idx = ismember(lower(updatedinfo.Subject),lower(keymouse)); % find index of mosekey sample name that matches data sample
    num_matches = sum(idx); % number of matching rows from data
    matching_data_row_num=find(idx); % row index numbers from data that have mouse
    condition = string(mousekey.Condition(x));
    % create a line here if the mouse has a different blinded name
      
    placeholder(matching_data_row_num) = condition;
    x = x+1;
end

updatedinfo.Condition=placeholder; % updated with condition data

%% check for outliers based on condition
temp=table;
c = 1;
allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
while (c<=length(allConditions))
    group=allConditions(c);
    select_idx = updatedinfo.Condition==group ; % logical indexing 
    selection = updatedinfo(select_idx,:);
    selection.(3) = filloutliers(selection.(3),NaN); %replace outliers from median with NaNs
    temp=[temp; selection];
    c=c+1;
end
updatedinfo=temp;
updated_unblindedinfo=temp;
% updated_unblindedinfo=updatedinfo;
updatedfile=([logpathname 'Data\' region '_' stringToBeFound '__appended_unblinded' date '.mat']);
save(updatedfile,'updated_unblindedinfo');
disp(['saved unblinded data ' updatedfile]);
save([logpathname 'Data\updated_unblindedinfo.mat'],'updatedinfo');





%% run graphing after unblinding above

clf; %clears current figure
disp('making graphs now');

%% get paths
addpath(genpath(pwd)); %code location
datapath = ([pwd '\Data\']);
addpath(datapath);
info=load('updated_unblindedinfo.mat');
updated_unblindedinfo=info.updatedinfo;

% set path to save to
figpathname = ([pwd '\Graphs\']);
addpath(figpathname);

% add violin plot folder to path
% addpath('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Violinplot-Matlab-master');
% addpath('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Violinplot-Matlab-master');
% addpath('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Violinplot-Matlab-master');
violinplot_location=('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_code\Violinplot-Matlab-master');

% "% plot the Violin", somewhere around line 145, you can add 'Linestyle', 'none' to the "fill" function and that will remove the outline
%% group Violin Plots
% figure;
updated_unblindedinfo.(3)=double(updated_unblindedinfo.(3));
%options ViolinColor, EdgeColor, MedianColor, BoxColor..
figure('Position', [19 402 462 398]);
vplot = violinplot(updated_unblindedinfo.(3), updated_unblindedinfo.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'});
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
ylim([0,(max(updated_unblindedinfo.(3))+10)]);
saveas(gcf,[figpathname stringToBeFound '_violin_plot_' date '.png']);

%% individual Violin Plots
% clf;
% figure;
% allConditions = {'Vehicle + 40 Hz', 'Vehicle + Light', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'};
% conditionColors = [0.9290 0.6940 0.1250 ; 0.6350 0.0780 0.1840 ; 0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330]; 
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
% vplot = violinplot(updated_unblindedinfo.(3), updated_unblindedinfo.Subject, 'BoxColor', [ 0 0 0]);
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
figure;
% Anova
[p,t,stats] = anovan(updated_unblindedinfo.(3), {updated_unblindedinfo.Condition}, "Varnames", "Condition");
manovafile = ([logpathname '\Stats\Anova_' stringToBeFound '_' date ]);
save(manovafile,'t');

% run posthoc comparisons between groups 
posthoc_comps = multcompare(stats, 'CType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
mmulitifile = ([logpathname '\Stats\TukeyKramerComps_' stringToBeFound '_' date]);
save (mmulitifile, 'posthoc_comps');
disp('Anova and post-hoc comparison stats saved');

%% Means and Stdev
% [G,conditions] = findgroups(updatedinfo.Condition);
% avgdata = splitapply(@mean,updatedinfo.(4),G);
% stdevdata = splitapply(@std,updatedinfo.(4),G);
% avgdata.(1) = stringToBeFound;

[means,sems] = grpstats(updated_unblindedinfo.(3),updated_unblindedinfo.Condition,["mean","sem"]);
avgdata=table(means, sems);
gnames={'Vehicle + 40Hz', 'Vehicle + Light', 'NFkBInhibitor + 40Hz', 'NFkBInhibitor + Light'};
avgdata.Properties.RowNames = gnames;
save(([logpathname '\Stats\Averages' stringToBeFound '_' date]), 'avgdata');
disp('averages and sems saved');