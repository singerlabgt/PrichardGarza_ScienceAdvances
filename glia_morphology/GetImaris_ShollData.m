%Get Imaris Sholl Data
clc;
clear;

%% set up paths run from Inhibitor Code Folder
 addpath(genpath(pwd)); %code location
%directory location for final data in Box
% logpathname = ('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
% logpathname = ('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
logpathname = ('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');

%location of data in Emma's Neurocloud Folder
datapathname = ('T:\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
% datapathname = ('Y:\singer\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
% datapathname = ('T:\singer\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
addpath(genpath(datapathname)); %data location including all folders and subfolders for each animal

% create a new info table to populate
originfo = array2table(NaN(1,5));
originfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'Radius' 'NumIntersections'};

%% edit below for region
region = 'VC';

%% **** specify string to find here **** run one at a time
stringToBeFound = 'No._Sholl_Intersections';
disp(['now calculating ' stringToBeFound ' for each sample in the ' region ' using means and maxes for each Filament ID...']);

%% read in folders and begin iteration through folder
filesAndfolders = dir(datapathname); % returns all files and folders in directory
numoffolders = length(filesAndfolders); % the number of folders in the directory
b = 1:size(filesAndfolders);
i=3; % set a folder to start, first 2 folders are filler '..' and should be skipped for this process

updatedinfo = originfo; %create table to update each iteration

while(i<=numoffolders)
    curdir = filesAndfolders(i).name;    % Store the name of the folder
    mouseID=string(curdir);
    disp(['finding morphology data for' mouseID]);
    mouseinfo = ([datapathname curdir]); %create a path to the folder
    mousedir = dir(mouseinfo); %make it a directory
    num_dir=numel(mousedir); %how many files and folders
%     all_filesandfoldersindir = dir(fullfile(mouseinfo,'**\**\*.*')); %count folders and files
    all_filesandfoldersindir = dir(fullfile(mouseinfo,'\**\*.*'));

    filesindir = all_filesandfoldersindir(~([all_filesandfoldersindir.isdir])); %get just files
    numberffiles = length(filesindir); % the number of files in the directory
    a=1; %create an interation
    while(a<=numberffiles) %start the while loop through files within the directory
        filename = filesindir(a).name;  % Store the name of the file
        found = strfind(filename, stringToBeFound); % look for specific variable in file name listed above
        if ~isempty(found) % if the results are not empty, it has found the file name
            foundString = strcat('Found in file ------', filename);
            disp(foundString); % sanity check
            thismouse = filename;
            data = readtable(thismouse);
%             data.(1) = filloutliers(data.(1),NaN); %replace outliers from median with NaNs
            num_cells = height(data);
            dummy_mouse = repmat(mouseID,num_cells,1);
            dummy_region = repmat(string(region),num_cells,1);
            newinfo = table(dummy_mouse, dummy_region, data.ID, data.Radius, data.FilamentNo_ShollIntersections);
            newinfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'Radius' 'NumIntersections' };
            updatedinfo = [updatedinfo; newinfo]; %append rows to table for each subset

            break;
        end
        a=a+1; % iterate through additional files in folder if not found
    end
    i = i+1; % move onto the next mouse ID
end        
%%        
disp(['updating info to table for ' stringToBeFound]);
disp(['in the region ' region]);
disp('done');

updatedinfo(1,:) = []; %delete first blank row
% updatedinfo=renamevars(updatedinfo,"VariableofInterest",stringToBeFound);
% updated_unblindedinfo=updatedinfo;

% pass to unblinding function but turn off graphing function at the end of
% unblinding function or will make weird violin plots..
unblinding_function(updatedinfo,region,stringToBeFound,logpathname)
load([logpathname 'Data\updated_unblindedinfo.mat']);
updated_unblindedinfo=updatedinfo;

%% check for outliers by radii
% updatedinfo.outs = filloutliers(updatedinfo.NumIntersections,NaN, 'movmedian',[3 3]); %replace outliers from median with NaNs
updated_unblindedinfo.Condition=categorical(updated_unblindedinfo.Condition);
conditionlist = unique(updated_unblindedinfo.Condition);
radiilist = unique(updated_unblindedinfo.Radius);
outliers_removed = updated_unblindedinfo(1,:);
k=1;
while k<=length(conditionlist)
    grp=conditionlist(k);
    r=0;
    sect=(updated_unblindedinfo.Condition == grp);
    sub_sect=updated_unblindedinfo(sect,:);
    while r<=height(radiilist)
        radi=(sub_sect.Radius == r);
        sub_sectradi=sub_sect(radi,:);
        numcells_at_radius=unique(sub_sectradi.Subject);
        if length(numcells_at_radius)<2
            r=r+1;
        else
            sub_sectradi.NumIntersections=filloutliers(sub_sectradi.NumIntersections,NaN, 'median');
            outliers_removed=[outliers_removed; sub_sectradi];
            r = r+1;
        end
    end
    k=k+1;
end

outliers_removed(1,:) = []; %delete first blank row
outliers_removed(outliers_removed.Radius>58,:)=[];

updatedfile=([logpathname 'Data\' region '_' stringToBeFound '__appended_' date '.mat']);
save(updatedfile,'updatedinfo');

% % pass to unblinding function but turn off graphing function at the end of
% % unblinding function or will make weird violin plots..
unblinding_function(updatedinfo,region,stringToBeFound,logpathname);

%% bin by 5
% % updatedinfo=updated_unblindedinfo;
% numradi=height(unique(updated_unblindedinfo.Radius));
% five_idx=0:5:75;
% five_idx';
% five_idx_rows=ismember(updated_unblindedinfo.Radius, five_idx);
% updatedinfo=updated_unblindedinfo(five_idx_rows,:);
%% find averages across radii
disp ('getting averages across radii');
meanbyradi = varfun(@mean,updatedinfo,'GroupingVariables',{'Radius', 'Condition' },'InputVariables',{'NumIntersections'});
%     meanbyradi = varfun(@mean,uptable,'GroupingVariables',{'Mouse', 'Radius', 'Condition' },'InputVariables',{'NumIntersections'});

nfile=([logpathname 'Data\' region '_Num_Sholl_Intersections_RadiiMeans']);
save(nfile,'meanbyradi');



%% Graph Data by Condition
clf; %clears current figure
disp(' ');
disp('making group graphs now');

% set path to save grpahs to
figpathname = ([pwd '\Graphs\']);
addpath(figpathname);

% same as  mean by radii above but reordered by condition first
[G,group,radi] = findgroups(meanbyradi.Condition,meanbyradi.Radius); %find groups
[group_intersectmeans] = splitapply(@mean,meanbyradi.mean_NumIntersections,G); %get averages for condition & radi

% !!!! check order of this below - does not match order of mean by radi or
% splitapply function
[means,sems] = grpstats(updatedinfo.NumIntersections,{updatedinfo.Radius, updatedinfo.Condition},["mean","sem"]);

result = table(group,radi,group_intersectmeans); %create a table of results
% 'GroupOrder', {'Vehicle + Light', 'Vehicle + 40 Hz', 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + Light'},

colors = [0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250]; %NFkB+40Hz, NFkB+Light, Vehicle+40Hz, Vehicle+Light
makegraph = gscatter(radi,group_intersectmeans,group, colors);
tempGroup=unique(group);
% hold on
% for igroup=1:length(tempGroup)
%     plotI=find(group==tempGroup(igroup));
%     plot(radi(plotI),group_intersectmeans(plotI),'-','color', colors(igroup,:));
% end
makegraph(4).MarkerFaceColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
makegraph(4).MarkerEdgeColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow
% makegraph(4).MarkerLineColor =[0.9290 0.6940 0.1250]; %vehicle + light = yellow

makegraph(3).MarkerFaceColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
makegraph(3).MarkerEdgeColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red
% makegraph(3).MarkerLineColor = [0.6350 0.0780 0.1840]; %vehicle + 40Hz = red

makegraph(1).MarkerFaceColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
makegraph(1).MarkerEdgeColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green
% makegraph(1).MarkerLineColor = [0.4660 0.6740 0.1880]; %NFkB inhibitor + 40Hz = green

makegraph(2).MarkerFaceColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
makegraph(2).MarkerEdgeColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue
% makegraph(2).MarkerLineColor = [0.3010 0.7450 0.9330]; % NFkB inhibitor + light lbue

set(makegraph, 'linestyle', '-');
xlabel('Sphere Radius (um)');
ylabel([region ' Average Number of Sholl Intersections']);
ylim([0,(max(group_intersectmeans)+1)]);
xlim([0, 60]);
saveas(gcf,[figpathname region '_AverageSholl_Intersections_scatterplot_' date '.png']);

%% set plot here to make graph w/o dots 
hold on
for igroup=1:length(tempGroup)
    plotI=find(group==tempGroup(igroup));
    plot(radi(plotI),group_intersectmeans(plotI),'-','color', colors(igroup,:), 'LineWidth', 2);
end
hold off
xlabel('Sphere Radius (um)');
ylabel([region ' Average Number of Sholl Intersections']);
legend({ 'NFkB Inhibitor + 40 Hz', 'NFkB Inhibitor + No Flicker', 'Vehicle + 40 Hz', 'Vehicle + No Flicker',})
saveas(gcf,[figpathname region '_AverageSholl_Intersections_plot_' date '.png']);

%% graph data by condition and by mouse (lots of lines!)
% clf; %clears current figure
disp(' ');
disp('making individual data point graphs now');

% get means by condition and by subject
meanbyradi = varfun(@mean,updatedinfo,'GroupingVariables',{'Subject', 'Radius', 'Condition' },'InputVariables',{'NumIntersections'});

[G,group,subject,radi] = findgroups(meanbyradi.Condition,meanbyradi.Subject, meanbyradi.Radius);
[group_intersectmeans] = splitapply(@mean,meanbyradi.mean_NumIntersections,G);
result = table(group,subject,radi,group_intersectmeans);
g = {result.group,result.subject};

figure;
makegraph = gscatter(radi,group_intersectmeans,g);
set(makegraph, 'linestyle', '-');
xlabel('Sphere Radius (um)');
ylabel([region ' Average Number of Sholl Intersections']);
ylim([0,(max(group_intersectmeans)+1)]);
saveas(gcf,[figpathname region '_AverageSholl_Intersections_byMouse_plot_' date '.png']);


