function imaging_time_tracesN=normalizeActivityByPreTrial(imaging_time_traces, t, st0, en0);
st=findClosestDouble(t,st0);
en=findClosestDouble(t,en0);

meanvals = mean(imaging_time_traces(:, st:en, :),2);
stdvals = std(imaging_time_traces(:, st:en, :),[],2);
imaging_time_tracesN = imaging_time_traces -repmat(meanvals, [1, size(imaging_time_traces,2), 1]);
imaging_time_tracesN = imaging_time_tracesN./repmat(stdvals, [1, size(imaging_time_tracesN,2), 1]);
