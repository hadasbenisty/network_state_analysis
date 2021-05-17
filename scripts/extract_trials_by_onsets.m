function [trials_data, trials_inds] = extract_trials_by_onsets(t_imaging, onset_timing, before_win, after_win, ...
    fsimaing, X)


if size(X, 1) == 1
    X = X.';
end
l = 1;trials_inds=[];
for onset_i = 1:length(onset_timing)
    if onset_timing(onset_i)<t_imaging(1)
        continue;
    end
    
    inds_imaging = findClosestDouble(t_imaging, onset_timing(onset_i));
    st_ind = inds_imaging - before_win*fsimaing;
    end_ind = inds_imaging + after_win*fsimaing;
    
    if st_ind < 1 
        continue;
    end
    if end_ind > size(X, 2)
        break;
    end
    trials_data(:, :, l) = X(:, st_ind:end_ind);
    trials_inds(end+1) = onset_i;
    l = l + 1;
end