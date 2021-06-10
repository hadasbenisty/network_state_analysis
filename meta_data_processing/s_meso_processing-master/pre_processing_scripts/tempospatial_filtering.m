function outsig = tempospatial_filtering(insig, signame, outputpathfilt, R, C)
addpath(genpath('../meta_data_processing/'))
addpath(genpath('../pre_processing_scripts/'));
addpath(genpath('../correlation/'));
filtlenDenoise=100;
wcutoff = 0.4;
nrnsBatchSz = 5e3;

patchsize = 2;
th = 0;


mkdir([outputpathfilt signame ]);  
disp(['Time filtering ' signame] );
outfilesBlueFilt = time_filteringBatch(insig, nrnsBatchSz, [outputpathfilt signame], filtlenDenoise, wcutoff);
disp(['Time filtering ' signame] );
outsig = nlmfilterbatch(outfilesBlueFilt, R, C, th, patchsize, [outputpathfilt signame]);
outsig=outsig(:,((filtlenDenoise/2)+1):end); %compensate for the delay in timing caused by fir filter 
end
function outblue = nlmfilterbatch(files, R, C, th, patchsize, signame)
for k = 1:length(files)
    tic;
    disp(['Loading batch no. ' num2str(k) ' of ' num2str(length(files))]);
    load(files{k}, 'data_time_filt','nrnsinds');
    blue_filt(nrnsinds, :) = data_time_filt; %#ok<AGROW>
    toc;
end
    outblue = local_non_localmean(1:R*C, R,C,blue_filt, patchsize, th, signame);
end
