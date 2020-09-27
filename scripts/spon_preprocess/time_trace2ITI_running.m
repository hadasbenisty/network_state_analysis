function [t_trials, trials_imaging_A] = time_trace2ITI_running(X, t_imaging, on_timestamps,fsimaing,after_win)
trials_imaging_A=[];t_trials=[];
for run_i = 1:length(on_timestamps)
    on_inds = findClosestDouble(t_imaging, on_timestamps(run_i));
    %off_inds = findClosestDouble(t_imaging, off_timestamps(run_i)); 
    if (on_inds+after_win*fsimaing) > size(X,2)        
    else
   temptrial=X(:, on_inds:on_inds+after_win*fsimaing);   
   temptimestamps=t_imaging(on_inds:on_inds+after_win*fsimaing);
   trials_imaging_A =cat(3, trials_imaging_A, temptrial);
   %trials_imaging_A(56,:,run_i)=temptrial;
   t_trials = cat(2, t_trials, temptimestamps);
    end
end
% third_stim = findClosestDouble(t_imaging, running_timestamps(3));
% plot(t_trials(1:third_stim).');
end

