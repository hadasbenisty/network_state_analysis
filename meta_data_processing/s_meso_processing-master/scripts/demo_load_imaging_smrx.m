addpath(genpath('../pre_processing_scripts'));
cedpath = '../pre_processing_scripts/utils/CEDS64ML';
addpath(genpath('../parcellation/'));
fsspike2 = 5e3;
smrxfile = '\\128.36.220.167\vivo\Lan\Meso-imaging\xs\xs_D1';
imagingfile = '\\128.36.220.173\vivo2\Lan\Meso-imaging\xs\meso_cleaned\xs_D1_PixxTime_dff.mat';
% load smrx
[timing, channels_data] = process_spike2(cedpath, '',smrxfile, fsspike2);

t_spike2 = linspace(0, length(channels_data.startsig)-1,length(channels_data.startsig))/fsspike2; 
% imaging
load(imagingfile);
% parcellation
[Allparcells, parcells] = getParcellsByLansAllansAtlas;

for par_i = 1:size(parcells.indicators,3)
    inds = find(parcells.indicators(:, :, par_i) == 1);
    dF_par(par_i, :) = mean(PixxTime_dff(inds, :));
end
% time line for imaging
t_imaging = timing.mesostart/fsspike2;
dF_par=dF_par(:, 1:length(t_imaging));
% sanity check to see that smrx and imaging data are alinged, plotting v1
% and stim traces. Zoom in to see that they match
figure;
h(1) = subplot(2,1,1);plot(t_imaging, dF_par(1,:));
xlabel('Time [sec]');ylabel('\Delta F/F');
title('V1 Response');
h(2) = subplot(2,1,2);plot(t_spike2, channels_data.startsig);
hold all;plot(t_spike2, channels_data.air_puff);
linkaxes(h,'x');
title('Smrx Data');legend('Stim','Air Puff');

stim_timestamps = timing.stimstart/fsspike2;
before_win = 120;
after_win = 139;
for stim_i = 1:length(stim_timestamps)
inds_imaging = findClosestDouble(t_imaging, stim_timestamps(stim_i));
imagingAl(:, :, stim_i) = dF_par(:, inds_imaging-before_win:inds_imaging+after_win);
end
fsimaing=33.3334;
t_trials = linspace(-before_win, after_win, after_win+before_win+1)/fsimaing;

figure;subplot(2,1,1);
imagesc(t_trials, 1:size(imagingAl,1), squeeze(imagingAl(1,:,:))')
xlabel('Time [sec]');ylabel('Trials');title('V1 Responses Vs Trials');
subplot(2,1,2);plot(t_trials, mean(imagingAl(1,:,:),3))
xlabel('Time [sec]');ylabel('\Delta F/F');
title('V1 Averaged Accross Trials');axis tight;



