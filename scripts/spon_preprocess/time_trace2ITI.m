function [t_trials, trials_imaging_A] = time_trace2ITI(X, t_imaging, stim_timestamps, before_win, after_win, fsimaing)
first_stim = findClosestDouble(t_imaging, stim_timestamps(1));
trials_imaging_A = X(:, 1:first_stim - before_win*fsimaing);
prev_stim = first_stim;
t_trials =t_imaging(1:first_stim - before_win*fsimaing);
for stim_i = 2:length(stim_timestamps)
    inds_imaging = findClosestDouble(t_imaging, stim_timestamps(stim_i));
    if (inds_imaging - before_win*fsimaing) > size(X,2)
    else
   trials_imaging_A =cat(2, trials_imaging_A, X(:, prev_stim+after_win*fsimaing:inds_imaging - before_win*fsimaing));
   t_trials = cat(1, t_trials, t_imaging(prev_stim+after_win*fsimaing:inds_imaging - before_win*fsimaing));
    end
 prev_stim =  inds_imaging;
end
% third_stim = findClosestDouble(t_imaging, running_timestamps(3));
% plot(t_trials(1:third_stim).');
end




