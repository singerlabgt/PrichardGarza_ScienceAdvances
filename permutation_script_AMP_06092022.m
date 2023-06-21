% permutations 
% read in data, change conditionVal, pass to permutation script,
% save out LV1 scores, compute means for each group on permuted scores
% compare to control group
clear all;
clear;
clc;

%% specify number of permutations
% R=1001; % Number of permutations
R = 1;
for permutecount = 1:R % set up iteration through number of permutations

    % read in data
    dataz=readtable('cytokineDataCorrected.xlsx');

    %dataPath='C:\Users\kgarza6\Box\Project_FlickerNeuroimmune_Team\Data\KristieRawFlickerData\Flicker 701-760_microgliadepeltion_Nov2020\cytokineDataCorrected.xlsx';

    %groupnames:
    GroupsforAnalysis=["40Hz";"40HzPLX";"20Hz"];

    %%Find Analysis Groups by comparing name of groups to group names in excel

    data1=dataz;
    data2=dataz;
    data3=dataz;
    data1(~ismember(dataz.conditionVal,GroupsforAnalysis(1)),:)=[]; % remove non group members
    data2(~ismember(dataz.conditionVal,GroupsforAnalysis(2)),:)=[]; % remove non group members
    data3(~ismember(dataz.conditionVal,GroupsforAnalysis(3)),:)=[]; % remove non group members

    data = [data1; data2; data3]; % combine values back into a data structure

    if permutecount == 1
        % run the original script once to get the true values first
        filename = 'cytokineDataCorrected.xlsx';
    else
    %% change labels of mice condition assignments to randomly permuted names
        permutated_new_indices = randperm(height(data))'; %to generate the indices you want your current real data to go to 
        permuted_data = data.conditionVal(permutated_new_indices); %to reorder
        data.conditionVal = permuted_data;
        filename = 'PLSDA_corrected_cytokine_data.xlsx';
        writetable(data,filename,'Sheet',1);
    end
    
% call PLSDA function here
    datamean = Fig3_permutations_PLSDA_manuscriptversion(filename); % feed in data with new labels, get out LV1means
    
    % return LV1 means called datamean rows are (20Hz, 40Hz, 40HzPLX)
    % column 1 = LV1 means, column 2 = LV2 means
  
    if permutecount == 1
        truemeans = array2table(datamean); % save out the origina, non permuted values
        truemeans.Properties.VariableNames = {'LV1', 'LV2'};
        truemeans.Properties.RowNames = {'20Hz', '40Hz', '40HzPLX'};
        save('truemeansL1L2.mat', 'truemeans');
        a = 'initial true means';
    elseif permutecount == 2
        means = array2table(datamean); % start creating a table of LV values
        a=num2str((permutecount-1));
        means.Properties.VariableNames = {['LV1 perm ' a], ['LV2 perm ' a]};
        means.Properties.RowNames = {'20Hz', '40Hz', '40HzPLX'};
    else
        newmeans = array2table(datamean);
        a=num2str((permutecount-1));
        newmeans.Properties.VariableNames = {['LV1 perm ' a], ['LV2 perm ' a]};
        means = [means newmeans]; % append new LV values to table each iteration
    disp(['permute count ' a]);
    end
end

disp (['completed ' a ' permutations']);
means.Properties.RowNames = {'20Hz', '40Hz', '40HzPLX'};

disp ('now check means table');
% Then compute means for each group on permuted scores and compare to control group
meansfile = 'PLSDA_permuted_cytokine_L1L2means.xlsx';
% writetable(means,meansfile,'Sheet',1); %turn off to not overwrite 1k permute file

truemeans= load('truemeansL1L2.mat','truemeans');

means = readtable(meansfile);
means.Properties.RowNames = {'20Hz', '40Hz', '40HzPLX'};

%% check permutations data on only odds for lv1, evens for lv2, otherwise will do both if on
%% lv2 = truemeans(1:end, 2:2:end) example
l2means = means(1:end, 2:2:end);
l2means_sqrd = table2array(l2means); %change to array
% l2means_sqrd = l2means_sqrd.^2; %square all values in array
array2table(l2means_sqrd);
l2truemeans = truemeans.truemeans(1:end, 2:2:end);

%% lv1 = truemeans(:,1:2:end) example
l1means = means(:,1:2:end);
l1means_sqrd = table2array(l1means); %change to array
% l1means_sqrd = l1means_sqrd.^2; %square all values in array
array2table(l1means_sqrd);
l1truemeans = truemeans.truemeans(:,1:2:end);

%% convert table values to array for math
l1_group20hz = table2array(l1means(1,:));
l1_group40hz = table2array(l1means(2,:));
l1_group40hzplx = table2array(l1means(3,:));

l2_group20hz = table2array(l2means(1,:));
l2_group40hz = table2array(l2means(2,:));
l2_group40hzplx = table2array(l2means(3,:));

%% find square root of differences between permuted means for groups for levels L1 and L2
% square root of (LV1_40Hz - LV1_20Hz)^2 + (LV2_40Hz - LV2_20Hz)^2
diff40vs20 = ((l1_group40hz - l1_group20hz).^2) + ((l2_group40hz - l2_group20hz).^2);
diff40vs20 = sqrt(diff40vs20);
% square root of (LV1_40Hz - LV1_40Hzplx)^2 + (LV2_40Hz - LV2_40Hzplx)^2
diff40vs40plx = ((l1_group40hz - l1_group40hzplx).^2) + ((l2_group40hz - l2_group40hzplx).^2);
diff40vs40plx = sqrt(diff40vs40plx);
% square root of (LV1_40Hzplx - LV1_20Hz)^2 + (LV2_40Hzplx - LV2_20Hz)^2
diff40plxvs20= ((l1_group40hzplx - l1_group20hz).^2) + ((l2_group40hzplx - l2_group20hz).^2);
diff40plxvs20 = sqrt(diff40plxvs20);

%% convert table values for true means to array for math
l1truemeans=table2array(l1truemeans);
l2truemeans=table2array(l2truemeans);

%% find true means differences
truediff40vs20 = ((l1truemeans(2,:) - l1truemeans(1,:))^2) + ((l2truemeans(2,:) - l2truemeans(1,:))^2);
truediff40vs20 = sqrt(truediff40vs20);

truediff40vs40plx = ((l1truemeans(2,:) - l1truemeans(3,:))^2) + ((l2truemeans(2,:) - l2truemeans(3,:))^2);
truediff40vs40plx = sqrt(truediff40vs40plx);

truediff40plxvs20 = ((l1truemeans(3,:) - l1truemeans(1,:))^2) + ((l2truemeans(3,:) - l2truemeans(1,:))^2);
truediff40plxvs20 = sqrt(truediff40plxvs20);

%% These lines count the number of times the difference in permuted means is greater than the true means and then averages to compute p values
pPerm1_40vs20 = mean(abs(diff40vs20) > abs(truediff40vs20)); % 40vs20 Hz
pPerm1_40vs20_N = sum(abs(diff40vs20) > abs(truediff40vs20)); % 40vs20 Hz

pPerm2_40vs40PLX = mean(abs(diff40vs40plx) > abs(truediff40vs40plx)); % 40vs40Hz PLX
pPerm2_40vs40PLX_N = sum(abs(diff40vs40plx) > abs(truediff40vs40plx)); % 40vs40Hz PLX

pPerm3_40PLXvs20 = mean(abs(diff40plxvs20) > abs(truediff40plxvs20)); % 40plxvs20
pPerm3_40PLXvs20_N = sum(abs(diff40plxvs20) > abs(truediff40plxvs20)); % 40plxvs20

%% combine and save file of pPerm values
pPerms = [pPerm1_40vs20 pPerm2_40vs40PLX pPerm3_40PLXvs20];
counts = [pPerm1_40vs20_N pPerm2_40vs40PLX_N pPerm3_40PLXvs20_N];
pPerms = [pPerms; counts];
pPerms = array2table(pPerms);
pPerms.Properties.VariableNames = {'40vs20', '40vs40PLX', '40PLXvs20'};
pPerms.Properties.RowNames = {'pPerm sig value', 'count of diff perms greater than true diff perms out of 1k perms'};
save('pPermsFinalValues.mat', 'pPerms');



%% Use this if only looking at one LV subtract the mean of the control group 20Hz from experimental groups (40Hz 40Hz PLX) 
% control_row_means = means(1:1,:); % first row is 20Hz
% control_row_perm_means = means("20Hz",:); % first row is labeled 20Hz
% forty_row_perm_means = means("40Hz",:); % second is labeled 40Hz
% forty_plx_row_perm_means = means("40HzPLX",:); % second is labeled 40Hz
% 
% meanOut1 = forty_row_perm_means{1:1,:} - control_row_perm_means{1:1,:}; %40hz-20hz_perm means
% meanOut1 = abs(meanOut1); % get absolute value
% % meanOut1 = array2table(meanOut1); %convert to table
% 
% meanOut2 = forty_plx_row_perm_means{1:1,:} - control_row_perm_means{1:1,:};%40hzplx-20hz_perm means
% meanOut2 = abs(meanOut2); % get absolute value
% % meanOut2 = array2table(meanOut2); %convert to table
% 
% %% subtract true means control group value from experimental true means groups
% controltrue = truemeans("20Hz",:);
% true40means = truemeans("40Hz",:);
% true40plxmeans = truemeans("40HzPLX",:);
% 
% meanTrue1 = true40means{1:1,:} - controltrue{1:1,:}; % 40hz-20hz_true means
% meanTrue1 = abs(meanTrue1); % get absolute value
% x = width(means); %find the numer of times to duplicate true means
% meanTrue1 = repmat(meanTrue1,1, x); % create a dummy true means table
% % meanTrue1 = array2table(meanTrue1);%convert to table
% 
% meanTrue2 = true40plxmeans{1:1,:} - controltrue{1:1,:}; %40hzplx-20hz_true means
% meanTrue2 = abs(meanTrue2); % get absolute value
% meanTrue2 = repmat(meanTrue2,1, x); % create a dummy true means table
% % meanTrue2 = array2table(meanTrue2);%convert to table

%% These lines count the number of times the difference in permuted means is greater than the true means and then averages to compute p values
% pPerm1 = mean(abs(meanOut1) > abs(meanTrue1)); % 40Hz returns
% pPerm1_N = sum(abs(meanOut1) > abs(meanTrue1));
% pPerm2 = mean(abs(meanOut2) > abs(meanTrue2)); % 40Hz PLX
% pPerm2_N = sum(abs(meanOut2) > abs(meanTrue2));
