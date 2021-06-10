%this code is to demix and align hannahs data for cohort 4
addpath(genpath('../network_state_analysis'))
csvfile = '../network_state_analysis/meta_data_processing/antara_exp/Processing_Pipeline_Full_expression_andHannah.xlsx';
processing_pipeline=readtable(csvfile);

%preprocess
processing_pipeline=processing_pipeline{strcmp(processing_pipeline{:,6},'Cohort 4'),1:7}
processing_pipeline{:,5};



for i=3:size(processing_pipeline,1)
    outputFolder=string(processing_pipeline{i,1});
    hannahnew_mecp21p_demixalign(outputFolder,char(processing_pipeline{i,2}),char(processing_pipeline{i,3}),char(processing_pipeline{i,4}),char(processing_pipeline{i,7}),char(processing_pipeline{i,6}))
    disp(num2str(i));
end