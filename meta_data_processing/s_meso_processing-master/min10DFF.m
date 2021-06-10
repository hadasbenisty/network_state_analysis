function [dff]=min10DFF(combinedDA) 
%get df/f where f0 is bottom 10% 
    dff=single(zeros(size(combinedDA,1), size(combinedDA,2)));
    parfor i =1:size(combinedDA,1) % for each pixel 
        subsample=combinedDA(i,:);     
        Ms = sort(subsample,'ascend');% Sort asending along time dimension
        F0Vals = Ms(1:ceil(length(Ms)*0.1)); % lower 10% of the values
        MeanF0=mean(F0Vals);
        tmp=((subsample-MeanF0)./MeanF0);
        dff(i,:)=tmp;
    end

