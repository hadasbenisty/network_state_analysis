function outfiles = time_filteringBatch(data, nrnsBatchSz, outputpath, filtlenDenoise, wcutoff)

b=fir1(filtlenDenoise,wcutoff);

nrnsBatchNum = floor(size(data, 1)/nrnsBatchSz);
outfiles = cell(nrnsBatchNum, 1);
for nrnsBatch_i = 1:nrnsBatchNum
    outfiles{nrnsBatch_i} = fullfile(outputpath, ['time_filt_data' num2str(nrnsBatch_i) '.mat']);
    if ~exist(outfiles{nrnsBatch_i}, 'file')
        disp(['Processing batch no. ' num2str(nrnsBatch_i) ' of ' num2str(nrnsBatchNum)]);
        tic;
        nrnsinds = 1 + (nrnsBatch_i - 1)*nrnsBatchSz:nrnsBatch_i*nrnsBatchSz;
        currBatch = double(data(nrnsinds, :));

        data_time_filt = filter(b,1,currBatch')';
        
        disp('saving');
        save(outfiles{nrnsBatch_i}, 'data_time_filt','nrnsinds','-V7.3');
        disp('Overall time per iteration: ');
        toc;
    end
end
if isempty(nrnsBatch_i)
    nrnsBatch_i=0;
end
nrnsBatch_i=nrnsBatch_i+1;
outfiles{nrnsBatch_i} = fullfile(outputpath, ['time_filt_data' num2str(nrnsBatch_i) '.mat']);
if ~exist(outfiles{nrnsBatch_i}, 'file')
    disp(['Processing batch no. ' num2str(nrnsBatch_i) ' of ' num2str(nrnsBatchNum+1)]);
    tic;
    nrnsinds = 1 + (nrnsBatch_i - 1)*nrnsBatchSz:size(data,1);
    currBatch = data(nrnsinds, :);


    data_time_filt = filter(b,1,currBatch')';
    disp('saving');
    save(outfiles{nrnsBatch_i}, 'data_time_filt','nrnsinds','-V7.3');
    disp('Overall time per iteration: ');
    toc;
end