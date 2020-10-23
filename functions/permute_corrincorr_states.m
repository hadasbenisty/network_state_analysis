function permute_corrincorr_states(x,i,animal)
%% correct vs incorrect, for  running

%concatenate conditions
concatenated_dff=cat(3,x.Correct_RunZscored,x.Incorrect_RunZscored);
%size of each condition and size of the entire dataset to permute
samplesize=size(x.Correct_RunZscored,3);
permuted_indices = randperm(size(concatenated_dff,3));
%save permutations
temp_cond1=concatenated_dff(:,:,permuted_indices(1:samplesize));
temp_cond2=concatenated_dff(:,:,permuted_indices(samplesize+1:length(permuted_indices)));

Correct_RunZscored=reshape(temp_cond1,size(temp_cond1,1),[]);
Incorrect_RunZscored=reshape(temp_cond2,size(temp_cond2,1),[]);

clearvars concatenated_dff samplesize permuted_indices temp_cond1 temp_cond2
%% correct vs incorrect, for NOT running

%concatenate conditions
concatenated_dff=cat(3,x.Correct_NotRunZscored,x.Incorrect_NotRunZscored);
%size of each condition and size of the entire dataset to permute
samplesize=size(x.Correct_NotRunZscored,3);
permuted_indices = randperm(size(concatenated_dff,3));
%save permutations
temp_cond1=concatenated_dff(:,:,permuted_indices(1:samplesize));
temp_cond2=concatenated_dff(:,:,permuted_indices(samplesize+1:length(permuted_indices)));

Correct_NotRunZscored=reshape(temp_cond1,size(temp_cond1,1),[]);
Incorrect_NotRunZscored=reshape(temp_cond2,size(temp_cond2,1),[]);
clearvars concatenated_dff samplesize permuted_indices temp_cond1 temp_cond2
%% correct vs incorrect, for low Pupil

%concatenate conditions
concatenated_dff=cat(3,x.Correct_LowPupZscored,x.Incorrect_LowPupZscored);
%size of each condition and size of the entire dataset to permute
samplesize=size(x.Correct_LowPupZscored,3);
permuted_indices = randperm(size(concatenated_dff,3));
%save permutations
temp_cond1=concatenated_dff(:,:,permuted_indices(1:samplesize));
temp_cond2=concatenated_dff(:,:,permuted_indices(samplesize+1:length(permuted_indices)));

Correct_LowPupZscored=reshape(temp_cond1,size(temp_cond1,1),[]);
Incorrect_LowPupZscored=reshape(temp_cond2,size(temp_cond2,1),[]);
clearvars concatenated_dff samplesize permuted_indices temp_cond1 temp_cond2
%% correct vs incorrect, for High Pupil

%concatenate conditions
concatenated_dff=cat(3,x.Correct_HighPupZscored,x.Incorrect_HighPupZscored);
%size of each condition and size of the entire dataset to permute
samplesize=size(x.Correct_HighPupZscored,3);
permuted_indices = randperm(size(concatenated_dff,3));
%save permutations
temp_cond1=concatenated_dff(:,:,permuted_indices(1:samplesize));
temp_cond2=concatenated_dff(:,:,permuted_indices(samplesize+1:length(permuted_indices)));

Correct_HighPupZscored=reshape(temp_cond1,size(temp_cond1,1),[]);
Incorrect_HighPupZscored=reshape(temp_cond2,size(temp_cond2,1),[]);
clearvars concatenated_dff samplesize permuted_indices temp_cond1 temp_cond2

save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',num2str(i),animal,'perm_state_trials_networkanalysis'),'Correct_RunZscored','Incorrect_RunZscored','Correct_NotRunZscored','Incorrect_NotRunZscored','Correct_LowPupZscored','Incorrect_LowPupZscored','Correct_HighPupZscored','Incorrect_HighPupZscored');

end

