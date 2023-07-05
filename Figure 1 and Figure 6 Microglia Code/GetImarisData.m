%GetImarisData
%Code to Retrieve and Colosolidate Morphology Data from Imaris CV files
clc;
clear;
clf;
%% set up paths
addpath(genpath(pwd)); %code location run from Inhibitor Code Folder
%directory location for final data in Box
% logpathname = ('C:\Users\Ashley\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
logpathname = ('C:\Users\aprichard3\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');
% logpathname = ('C:\Users\ashle\Box\Project_FlickerNeuroimmune_Team\Paper 2\Manuscript Data and Code\Inhibitor_Code\');

%location of data in Emma's Neurocloud Folder
% datapathname = ('T:\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
datapathname = ('Y:\singer\UndergradProjects\Emma\Microglia reconstuction\Post Processed Data Files\');
addpath(genpath(datapathname)); %data location including all folders and subfolders for each animal

% create a new info table to populate
originfo = array2table(NaN(1,4));
originfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'VariableofInterest'};
% use this line below if you have multiple samples for each mouse
% originfo = array2table(NaN(1,7));
% originfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};

%% edit below for region
region = 'VC';

%% morphology specifics to look for in the Imaris file name
% Area Measures (total cell area, soma area, average branch length) 
% Volume Measures (total cell volume, soma volume, average branch volume)
% microglia_number = 'Soma Area ImageJ.xlsx';
total_process_area = 'Filament_Area_(sum)'; % does not include soma
% soma_area = 'Soma Area Image J.xlsx';
branching_depth = 'Filament_Full_Branch_Depth';
total_process_length = 'Filament_Length_(sum)'; % length per process

total_process_volume = 'Filament_Volume_(sum)'; % total process volume per microglia - the soma (no soma here)
soma_volume = 'Soma_Volume';
total_process_volume = 'Filament_Dendrite_Volume_(sum)'; % total volume for each cell  not average
node_points = 'Filament_No._Dendrite_Branch_Pts'; % number of nodes
number_of_branches = 'Filament_No._Dendrite_Branches'; %number of branches resulting from nodes
average_segment_volume = 'Dendrite_Volume';
average_segment_length = 'Dendrite_Length';

%% **** specify string to find here **** run one at a time
stringToBeFound = average_segment_volume;
disp(['now calculating ' stringToBeFound ' for each sample in the ' region ' using means and maxes for each Filament ID...']);

%% read in folders and begin iteration through folder
filesAndfolders = dir(datapathname); % returns all files and folders in directory
numoffolders = length(filesAndfolders); % the number of folders in the directory
b = 1:size(filesAndfolders);
i=3; % set a folder to start, first 2 folders are filler '..' and should be skipped for this process

averages_table=table;
updatedinfo = originfo; %create table to update each iteration

while(i<=numoffolders)
    curdir = filesAndfolders(i).name;    % Store the name of the folder
%     mouseID = extractBetween(curdir,"(",")"); % get animal ID
%     mouseID=string(mouseID);
    mouseID=string(curdir);
    disp(['finding morphology data for' mouseID]);
%     if contains(region,'VC') % only need if multiple samples per animal
%         sampleID= extractBetween(curdir,"40x_","CA1"); %get specimen number
%         sampleID=string(sampleID);
%         sampleID= extractBetween(sampleID,"_","_"); 
%     end
    mouseinfo = ([datapathname curdir]); %create a path to the folder
    mousedir = dir(mouseinfo); %make it a directory
    num_dir=numel(mousedir);
%     filesindir = mousedir(~([mousedir.isdir])); % returns only files in the directory
%     all_filesandfoldersindir = dir(fullfile(mouseinfo,'**\**\*.*'));
    all_filesandfoldersindir = dir(fullfile(mouseinfo,'\**\*.*'));

    filesindir = all_filesandfoldersindir(~([all_filesandfoldersindir.isdir]));
%     a = subdir(fullfile(matlabroot, 'toolbox', 'matlab', '*.mat'))
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
            % find if there are outliers within the individual mouse data
%             data.(1) = filloutliers(data.(1),NaN); %replace outliers from median with NaNs
            
            num_cells = height(data);
            dummy_mouse = repmat(mouseID,num_cells,1);
            dummy_region = repmat(string(region),num_cells,1);
            newinfo = table(dummy_mouse, dummy_region, data.ID, data.(1));
            newinfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'VariableofInterest' };
%% use below to get means by groups for dendrite process length
            avgcalc= strcmp(stringToBeFound,'Dendrite_Length');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Filament_maxlens = splitapply(@max,data.DendriteLength,G); % find maximum dendrite length for each filament ID
                
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'DendriteLength'});
                meanbycell.MaxDendriteLength = Filament_maxlens;
                h = height(meanbycell);
                temp = {mouseID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};
                averages_table=[averages_table; avginfo];
                updated_avgfile=([logpathname 'Data\' region '_averages_' stringToBeFound '_appended_' date '.mat']);
                save(updated_avgfile,'averages_table');
                disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data\' stringToBeFound '_averages_table.mat'],'averages_table');
            end
%% get average data for dendrite volume
            avgcalc= strcmp(stringToBeFound,'Dendrite_Volume');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Filament_maxvols = splitapply(@max,data.DendriteVolume,G); % find maximum dendrite length for each filament ID
                
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'DendriteVolume'});
                meanbycell.MaxDendriteVolume = Filament_maxvols;
                h = height(meanbycell);
                temp = {mouseID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessVolume' 'MaxProcessVolume'};
                averages_table=[averages_table; avginfo];

                updated_avgfile=([logpathname 'Data\' region '_averages_' stringToBeFound '_appended_' date '.mat']);
                save(updated_avgfile,'averages_table');
                disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data\' stringToBeFound '_averages_table.mat'],'averages_table');
            end
%% use below to get means by groups for total process length
%             avgcalc= strcmp(stringToBeFound,'Filament_Length_(sum)');
%             if avgcalc == 1
%                 [G,maxlen] = findgroups(data.FilamentID);
%                 Filament_maxlens = splitapply(@max,data.FilamentLength_sum_,G); % find maximum dendrite length for each filament ID
%                 
%                 omean = @(x) mean (x,'omitnan');
%                 meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'FilamentLength_sum_'});
%                 meanbycell.MaxDendriteLength = Filament_maxlens;
%                 h = height(meanbycell);
%                 temp = {mouseID, region}; % create a table of new info
%                 avginfo = repmat(temp, h,1);
%                 avginfo = cell2table(avginfo);
%                 avginfo = [avginfo meanbycell];
%     % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
%                 avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};
%                 averages_table=[averages_table; avginfo];
%                 updated_avgfile=([logpathname 'Data\' region '_averages_' stringToBeFound '_appended_' date '.mat']);
%                 save(updated_avgfile,'averages_table');
%                 disp(['saved averages data ' updated_avgfile]);
%                 save([logpathname 'Data\' stringToBeFound '_averages_table.mat'],'averages_table');
%             end


            updatedinfo = [updatedinfo; newinfo]; %append rows to table for each subset
            
%             curmouse = (mouseID + sampleID); %get current sample name
%             disp (['sample' curmouse]);
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
updatedinfo=renamevars(updatedinfo,"VariableofInterest",stringToBeFound);
% updatedinfo.Filament_Full_Branch_Depth=updatedinfo.VariableofInterest;

updatedfile=([logpathname 'Data\' region '_' stringToBeFound '__appended_' date '.mat']);
save(updatedfile,'updatedinfo');

save([logpathname 'Data\updatedinfo.mat'],'updatedinfo');
disp(['saved ' updatedfile]);
disp('passing to unblinding function');  

unblinding_function(updatedinfo,region,stringToBeFound,logpathname)
