function corrmat_sh = extract_shuffled_segs(simname, segments_arousals, t_imaging, parcels_traces, state1, state2, REPS)

y1 = extract_segment(t_imaging, parcels_traces, segments_arousals.(state1));
y2 = extract_segment(t_imaging, parcels_traces, segments_arousals.(state2));
Y = [y1 y2];
labels = [ones(size(y1,2),1) ; 2*ones(size(y2,2),1)];
for r = 1:REPS
    inds = randperm(length(labels));
    labelsR = labels(inds);
    corrmat_sh.(state1)(:,:,r) = measure_weights_bysegs(1:sum(labelsR==1),Y(:, labelsR==1), simname);
    corrmat_sh.(state2)(:,:,r) = measure_weights_bysegs(1:sum(labelsR==2),Y(:, labelsR==2), simname);
end
