function plot_corr_by_state(animals, state1, state2, h, name1, name2)
parcels_names = get_allen_meta_parcels;
clr = colormap(redblue);
clr=cat(1, [.5 .5 .5], clr);

figure;
for ci=1:length(animals.type_lut)
subplot(length(animals.type_lut),3,ci);plot_corr_mat(nanmean(state1(:,:,animals.type_list==ci),3),parcels_names,[0 1]);
title([animals.type_lut{ci} ' ' name1]);
subplot(length(animals.type_lut),3,ci+3);plot_corr_mat(nanmean(state2(:,:,animals.type_list==ci),3),parcels_names,[0 1]);
title([animals.type_lut{ci} ' ' name2]);
A = nanmean(state1(:,:,animals.type_list==ci),3)-nanmean(state2(:,:,animals.type_list==ci),3);
A(h(:,:,ci)==1) = nan;
o=tril(ones(size(A)));
A(o==0) = 0;
subplot(length(animals.type_lut),3,ci+6);plot_corr_mat((A),parcels_names,[-1 1]*0.1);
title([animals.type_lut{ci} ' ' name1 '-' name2]);
end
colormap(clr);