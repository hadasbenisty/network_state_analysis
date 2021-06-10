function outfiles = saveBatchedNrnsByTime(data, timeBatchSz, filename, outputpath)
timeBatchNum = floor(size(data, 2)/timeBatchSz);
outfiles = cell(timeBatchNum, 1);
MAX_TRY = 100;
for timeBatch_i = 1:timeBatchNum
    outfiles{timeBatch_i} = fullfile(outputpath, [filename num2str(timeBatch_i) 'time.mat']);
    if ~exist(outfiles{timeBatch_i}, 'file')
        disp(['Saving batch no. ' num2str(timeBatch_i) ' of ' num2str(timeBatchNum)]);
        tic;
        timeinds = 1 + (timeBatch_i - 1)*timeBatchSz:timeBatch_i*timeBatchSz;
        currBatch = data(:, timeinds);
        
        
        
        disp('saving');
        suc = per_save(outfiles{timeBatch_i}, currBatch, timeinds);
        disp('Overall time per iteration: ');
        toc;
        if suc==0
            disp('Major problem, could not save');
            error('Major problem, could not save');
        end
        
    end
end
timeBatch_i=timeBatch_i+1;
timeinds = 1 + (timeBatch_i - 1)*timeBatchSz:size(data,2);
if isempty(timeinds)
    return;
end
outfiles{timeBatch_i} = fullfile(outputpath, [filename num2str(timeBatch_i) 'time.mat']);
if ~exist(outfiles{timeBatch_i}, 'file')
    disp(['Saving batch no. ' num2str(timeBatch_i) ' of ' num2str(timeBatchNum+1)]);
    tic;
    
    currBatch = data(:, timeinds);
    
    
    disp('saving');
    suc = per_save(outfiles{timeBatch_i}, currBatch, timeinds);
    if suc==0
        disp('Major problem, could not save');
        error('Major problem, could not save');
    end
    disp('Overall time per iteration: ');
    toc;
end