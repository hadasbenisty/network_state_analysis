function [pupil_t, pupil_trace] = time_trace2trials(pupil_Norm, pupil_time, t_imaging, stim_timestamps,  fspupilcam)
    bef = -t_imaging(1);

dur = round((bef+t_imaging(end))*fspupilcam);
    for k=1:length(stim_timestamps)
       tind = findClosestDouble(pupil_time, stim_timestamps(k)-bef); 
       pupil_trace(1,:,k) = pupil_Norm(tind:tind+dur);
    end
    pupil_t = linspace(-bef, dur/fspupilcam-bef, size(pupil_trace,2));

    
% first_stim = findClosestDouble(t_imaging, stim_timestamps(1));
% trials_imaging_A = X(:, 1:first_stim - before_win*fsimaing);
% prev_stim = first_stim;
% t_trials =t_imaging( 1:first_stim - before_win*fsimaing);
% for stim_i = 2:length(stim_timestamps)
%     inds_imaging = findClosestDouble(t_imaging, stim_timestamps(stim_i));
%    trials_imaging_A =...
%  cat(2, trials_imaging_A, X(:, prev_stim+after_win*fsimaing:inds_imaging - before_win*fsimaing));
% t_trials = cat(2, t_trials, t_imaging(prev_stim+after_win*fsimaing:inds_imaging - before_win*fsimaing));
% 
%    prev_stim =  inds_imaging;
% end

