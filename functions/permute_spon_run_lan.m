function permute_spon_run_lan(x,i,animal)
%concatenate conditions
concatenated_dff=cat(3,x.locomotionperiods,x.quiescenceperiods);
concatenated_timestamps=cat(2,x.t_wheelon,x.t_wheeloff);
%size of each condition and size of the entire dataset to permute
samplesize=size(x.locomotionperiods,3);
permutation = randperm(size(concatenated_dff,3));
%save permutations
temp_runningalldata = concatenated_dff(:,:,permutation(1:samplesize));
temp_notrunningalldata=concatenated_dff(:,:,permutation(samplesize+1:length(permutation)));
temp_runningalldata_t=concatenated_timestamps(:,permutation(1:samplesize));
temp_notrunningalldata_t=concatenated_timestamps(:,permutation(samplesize+1:length(permutation)));

runningalldata=reshape(temp_runningalldata,size(temp_runningalldata,1),[]);
notrunningalldata=reshape(temp_notrunningalldata,size(temp_notrunningalldata,1),[]);
runningalldata_t=reshape(temp_runningalldata_t,size(temp_runningalldata_t,1)*size(temp_runningalldata_t,2),1);
notrunningalldata_t=reshape(temp_notrunningalldata_t,size(temp_notrunningalldata_t,1)*size(temp_notrunningalldata_t,2),1);

save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',num2str(i),animal,'perm_running_time_traces'),'runningalldata','runningalldata_t','notrunningalldata_t','notrunningalldata');
end



