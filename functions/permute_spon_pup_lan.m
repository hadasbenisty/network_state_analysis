function permute_spon_pup_lan(x,i,animal)
%concatenate conditions
concatenated_dff=cat(3,x.puphigh_on,x.puplow_on);
concatenated_timestamps=cat(2,x.t_puphigh_on,x.t_puplow_on);
%size of each condition and size of the entire dataset to permute
samplesize=size(x.puphigh_on,3);
permutation = randperm(size(concatenated_dff,3));
%save permutations
temp_highpupalldata = concatenated_dff(:,:,permutation(1:samplesize));
temp_lowpupalldata=concatenated_dff(:,:,permutation(samplesize+1:length(permutation)));
temp_highpupalldata_t=concatenated_timestamps(:,permutation(1:samplesize));
temp_lowpupalldata_t=concatenated_timestamps(:,permutation(samplesize+1:length(permutation)));

pupilhighalldata=reshape(temp_highpupalldata,size(temp_highpupalldata,1),[]);
pupillowalldata=reshape(temp_lowpupalldata,size(temp_lowpupalldata,1),[]);
pupilhighalldata_t=reshape(temp_highpupalldata_t,size(temp_highpupalldata_t,1)*size(temp_highpupalldata_t,2),1);
pupillowalldata_t=reshape(temp_lowpupalldata_t,size(temp_lowpupalldata_t,1)*size(temp_lowpupalldata_t,2),1);

save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',num2str(i),animal,'perm_pupilhigh_pupillow'),'pupilhighalldata','pupilhighalldata_t','pupillowalldata_t','pupillowalldata');
end

