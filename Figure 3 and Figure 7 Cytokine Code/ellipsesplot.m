%create Ellipses
function= ellipsesplot(BarLabels,treatment,region.pls.T)
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
    
     %create covariance matrix
    group1covmatrix=cov(GroupLVscores{1,1});
    group2covmatrix=cov(GroupLVscores{1,2});
    group3covmatrix=cov(GroupLVscores{1,3});
    group4covmatrix=cov(GroupLVscores{1,4});
   
    %plot ellipses

error_ellipse(group1covmatrix,datamean(1,:),.95,'style','b')
error_ellipse(group2covmatrix,datamean(2,:),.95,'style','r')
error_ellipse(group3covmatrix,datamean(3,:),.95,'style','k')
error_ellipse(group4covmatrix,datamean(4,:),.95)
