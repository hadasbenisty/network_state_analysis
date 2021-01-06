function main_time_traces_to_trials




addpath(genpath('../meta_data_processing/'));
% load
[~, allen_parcels, maskAllen, maskAllenFrontal] = getParcellsByLansAllansAtlas;


% fltstr = 'pixelwise';
fltstr = 'spt';
%
animalNames = {'xt' 'xs' 'xu'  'xx'  'xz' 'xw'    };% 'xz' 'xx'
for ai = 1:length(animalNames)
    animal = animalNames{ai};
    [~,days] = animaltodays(animal);
    fsspike2=5e3;
    
    spike2pth0 = 'D:\datasets\meso\';
    
    dbstop if error;
    
    
    data_smr_path = 'X:\Lan\Meso-imaging\';
    cedpath = '../pre_processing_scripts/utils/CEDS64ML';
    
        switch animal
            case {'xu','xv','xt','xs'}
                fsimaing=33;
                delay_filt = 500;
            otherwise
                fsimaing=10;
                delay_filt=150;
        end
        for day_i = 1:length(days)
            disp([animal  ' '  num2str(days(day_i))]);
            
            datapath = ['X:\Hadas\Meso-imaging\lan\' animal 'psych\' fltstr '\'];
            spike2pth = fullfile(spike2pth0, animal);
         
            if exist(fullfile(datapath, [animal num2str(days(day_i)) 'imaging_time_traces_global_dfff.mat']), 'file')
                continue;
            end
            % spike2
            if ~exist(fullfile(spike2pth, ['spike2data',animal num2str(days(day_i)) '.mat']),'file')
                
                filename = fullfile(data_smr_path, animal, [animal '_D' num2str(days(day_i))]);
                if ~exist([filename '.smrx'],'file')
                    disp('No smrx');
                    continue;
                end
                [timing, channels_data] = process_spike2(cedpath, '',filename, fsspike2);
                t_spike2 = linspace(0, length(channels_data.startsig)-1,length(channels_data.startsig))/fsspike2;
                t_imaging = timing.mesostart/fsspike2;
                
                save(fullfile(spike2pth, ['spike2data',animal num2str(days(day_i)) '.mat']), 't_spike2', 'channels_data', 't_imaging',...
                    'timing', 't_spike2');
            else
                load(fullfile(spike2pth, ['spike2data',animal num2str(days(day_i)) '.mat']),'channels_data',...
                    'timing', 't_imaging');
                
            end
            %         t_imaging = timing.mesostart/fsspike2;
            before_win = 120;
            after_win = 139;
            % load parcells
            
            parfile = fullfile(datapath, [animal '_' num2str(days(day_i)) '_allen_dfff.mat']);
            parfileGal = fullfile(datapath,'gal',[animal '_' num2str(days(day_i)) ...
                '_global_dfff.mat']);
            parfileGrid = fullfile(datapath,[animal '_' num2str(days(day_i)) ...
                '_grid4_dfff.mat']);
            
            if ~exist(parfile, 'file')||~exist(parfileGal, 'file')|| ~exist(parfileGrid, 'file')       
                disp('files missing');
                continue;
            end
            pardataAllan = load(parfile);
            pardata_gal = load(parfileGal);
            
            [parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
            regionLabel.Allen = allen_parcels.regionNum;
            regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
            % [roiLabelsbyAllen_gal, regionLabel.Gal, maskByAllen_gal, maskByAllen.Gal] = get_gal_parcels_lables(animal);
            [parcels_names.Gal, finalindex.Gal, maskByAllen.Gal, regionLabel.Gal, roiLabelsbyAllen.Gal] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
                regionLabel.Allen, animal);
            pardata_grid = load(parfileGrid);
            [parcels_names.grid4, regionLabel.grid4,finalindex.grid4, ~, maskByAllen.grid4, roiLabelsbyAllen.grid4] = getAllenClusteringLabelsGrid(pardata_grid.par_inds, 4);
            roiLabelsbyAllen.Allen = 1:length(allen_parcels.names);
            
            
           
            
            %         delay_filt = 50;
            if fsimaing < 33
                before_win = 40;
                after_win = 45;
                t_imaging=t_imaging(1:2:end);
            end
            t_imaging = t_imaging(1:end-delay_filt);
            if length(t_imaging) > size(pardataAllan.parcels_time_trace,2)
                t_imaging=t_imaging(1:size(pardataAllan.parcels_time_trace,2));
            end
            X.Allen =pardataAllan.parcels_time_trace(:, 1:length(t_imaging));
            X.Gal =pardata_gal.parcels_time_trace(:, 1:length(t_imaging));
            X.grid4 =pardata_grid.parcels_time_trace(:, 1:length(t_imaging));
            names = fieldnames(X);
            for ni=1:length(names)
                X.(names{ni}) = X.(names{ni})(finalindex.(names{ni}),:);                
            end
            
            stim_timestamps = timing.stimstart/fsspike2;
            stim_timestamps=stim_timestamps(1:75);
            
            for stim_i = 1:length(stim_timestamps)
                tind = findClosestDouble(t_imaging, stim_timestamps(stim_i));
                if any(tind-before_win:tind+after_win > length(t_imaging))
                    for ni=1:length(names)
                imaging_time_traces.(names{ni})(:, :, stim_i) = nan(size(X.(names{ni}),1), length(tind-before_win:tind+after_win));
                    end
                else
                for ni=1:length(names)
                imaging_time_traces.(names{ni})(:, :, stim_i) = X.(names{ni})(:, tind-before_win:tind+after_win);
                imaging_time_traces.t = t_imaging(tind-before_win:tind+after_win)-stim_timestamps(stim_i);

                end
                end
            end
            
            if size(imaging_time_traces.Allen,3) ~= size(imaging_time_traces.Gal,3)
                error('check this');
            end
            trialsinds = 1:size(imaging_time_traces.Gal,3);
            blinkfile = fullfile('X:\Lan\Meso-imaging\', animal, [animal '_D' ...
                num2str(days(day_i))  '_blinksummary.mat']);
            if exist(blinkfile, 'file')
                trialslabels = load(blinkfile, 'blinksummary');
            else
                error('no blink file');
            end
            
            
            %% sti file
            stifile = fullfile('X:\Lan\Meso-imaging\', animal, [animal '_D' ...
                num2str(days(day_i))  '_sti.mat']);
            if ~exist(stifile, 'file')
                disp('no sti');
                continue;
            end
            sti = load(stifile,  'stiparameter');
            
            trialslabels.injectioncond = sti.stiparameter(trialsinds,end);
            trialslabels.contrastLabels = sti.stiparameter(trialsinds,1);
            trialslabels.blinksummary = trialslabels.blinksummary(trialsinds,1);
           
        
            figure;subplot(2,2,1);imagesc(imaging_time_traces.t, 1:75, squeeze(imaging_time_traces.Allen(1,:,trialslabels.blinksummary<3))');
            subplot(2,2,3);plot(imaging_time_traces.t, mean(imaging_time_traces.Allen(1,:,trialslabels.blinksummary<3),3));
            
            
            subplot(2,2,2);imagesc(imaging_time_traces.t, 1:75, squeeze(mean(imaging_time_traces.Gal(roiLabelsbyAllen.Gal==1,:,trialslabels.blinksummary<3),1))');
            subplot(2,2,4);plot(imaging_time_traces.t, mean(mean(imaging_time_traces.Gal(roiLabelsbyAllen.Gal==1,:,trialslabels.blinksummary<3),3),1));
            saveas(gcf, fullfile(datapath,[animal num2str(days(day_i))  '_V1timetraces_global_dfff.jpg']));
            close(gcf);
            save(fullfile(datapath, [animal num2str(days(day_i)) 'imaging_time_traces_global_dfff.mat']), 'imaging_time_traces' ,...
                'trialslabels', 'roiLabelsbyAllen', 'maskByAllen', 'regionLabel');
            clear imaging_time_traces;
            clear parcelsLabels;
            clear trialslabels;
            if exist(fullfile('X:\Hadas\Meso-imaging\lan\', [animal 'psych'], 'spike2Features', [animal num2str(days(day_i)) 'spike2.mat']),'file')
                continue;
            end
            t_spike2 = [1:length(channels_data.wheel)];
            [ond,ofd] = squaredetect(channels_data.diode,.5);
            
            for stim_i = 1:length(stim_timestamps)
                tind = findClosestDouble(t_spike2/fsspike2, stim_timestamps(stim_i));
                spike2time_traces.diode(:, :, stim_i) = channels_data.diode(tind-round(before_win/fsimaing*fsspike2):tind+round(after_win/fsimaing*fsspike2));
                spike2time_traces.air_puff(:, :, stim_i) = channels_data.air_puff(tind-round(before_win/fsimaing*fsspike2):tind+round(after_win/fsimaing*fsspike2));
                spike2time_traces.wheelspeed(:, :, stim_i) = channels_data.wheelspeed(tind-round(before_win/fsimaing*fsspike2):tind+round(after_win/fsimaing*fsspike2));
                spike2time_traces.wheel(:, :, stim_i) = channels_data.wheel(tind-round(before_win/fsimaing*fsspike2):tind+round(after_win/fsimaing*fsspike2));
                spike2time_traces.t(:, stim_i) = t_spike2(tind-round(before_win/fsimaing*fsspike2):tind+round(after_win/fsimaing*fsspike2));
            end
            

            mkdir(fullfile('X:\Hadas\Meso-imaging\lan\', [animal 'psych'], 'spike2Features'));
            save(fullfile('X:\Hadas\Meso-imaging\lan\', [animal 'psych'], 'spike2Features', [animal num2str(days(day_i)) 'spike2.mat']),  'spike2time_traces');
            clear spike2time_traces;
            
            
            %     outputfigs = 'C:\Users\Hadas Ben Esti\Dropbox (HigleyLab)\HigleyLab Team Folder\Hadas\meso\behave';
            
            %         [movData, tmov] = loadBehaveVideo(['X:\Lan\Meso-imaging\' animal '\' animal '_D' num2str(num2str(days(day_i))) '_*.mp4'],  before_win, after_win, stim_timestamps, t_behave_mov, outputfigs);
            
            % facemap
            %         facemapfile = dir(['X:\Hadas\Meso-imaging\lan\' animal '\videos\D' num2str(days(day_i)) '\*.mat']);
            %         if isempty(facemapfile)
            %             disp('no facemap');
            %         else
            %             load(fullfile(facemapfile.folder, facemapfile.name));
            %             facemap.pupil.area = proc.pupil.area;
            %             facemap.pupil.area_raw = proc.pupil.area_raw;
            %             facemap.svd = proc.motSVD{2};
            %               channels_data.pupil_frame=(channels_data.pupil_frame-nanmin(channels_data.pupil_frame)>0.5);
            %     [pupon,pupoff]=squaredetect(channels_data.pupil_frame,0.05);
            %     pupon=pupon(1:length(pupoff));
            %     t_behave_mov = pupon/fsspike2;
            %             [facemap_time_trace.t, facemap_time_trace.svd] = time_trace2trials(facemap.svd', t_behave_mov, stim_timestamps, before_win, after_win, fsimaing);
            %             [~, facemap_time_trace.pupil] = time_trace2trials(facemap.pupil.area', t_behave_mov, stim_timestamps, before_win, after_win, fsimaing);
            % save(fullfile('X:\Hadas\Meso-imaging\lan\',animal, 'videoFeatures', [animal num2str(days(day_i)) 'facemap.mat']), ...
            %             'facemap_time_trace');
            %         end
            % whisking
            wskile = fullfile('X:\Lan\Meso-imaging\', animal, [animal '_D' ...
                num2str(days(day_i))  '_whisker_clean.mat']);
            if ~exist(wskile, 'file')
                disp('no wskile');
                continue;
            end
            wsk = load(wskile,  'whisk_binary');
            wsk.whisk_binary =wsk.whisk_binary(1:min(length(wsk.whisk_binary),length(t_imaging)));
            for stim_i = 1:length(stim_timestamps)
                tind = findClosestDouble(t_imaging, stim_timestamps(stim_i));
                if any(tind-before_win:tind+after_win > length(wsk.whisk_binary))
                   
               wsk_time_trace.binary(:, :, stim_i) = nan(1, length(tind-before_win:tind+after_win));
                  
                else
                wsk_time_trace.binary(:, :, stim_i) = wsk.whisk_binary(tind-before_win:tind+after_win);
%                 imaging_time_traces.t = t_imaging(tind-before_win:tind+after_win)-stim_timestamps(stim_i);
                end
            end
            
            wskile = fullfile('X:\Lan\Meso-imaging\', animal, [animal '_D' ...
                num2str(days(day_i))  '_whisker.mat']);
            if ~exist(wskile, 'file')
                disp('no wskile');
                continue;
            end
            wsk = load(wskile,  'whisker');
            wsk.whisker =wsk.whisker(1:min(length(wsk.whisker),length(t_imaging)));
            
%             [wsk_time_trace.t, wsk_time_trace.whisker] = time_trace2trials(wsk.whisker', t_imaging, stim_timestamps(trialsinds), before_win, after_win, fsimaing);
            for stim_i = 1:length(stim_timestamps)
                tind = findClosestDouble(t_imaging, stim_timestamps(stim_i));
                if any(tind-before_win:tind+after_win > length(wsk.whisker))
                   
               wsk_time_trace.whisker(:, :, stim_i) = nan(1, length(tind-before_win:tind+after_win));
                  
                else
                wsk_time_trace.whisker(:, :, stim_i) = wsk.whisker(tind-before_win:tind+after_win);
                wsk_time_trace.t = t_imaging(tind-before_win:tind+after_win)-stim_timestamps(stim_i);
                end
            end
            mkNewDir(fullfile('X:\Hadas\Meso-imaging\lan\',[animal 'psych'], 'videoFeatures'))
            save(fullfile('X:\Hadas\Meso-imaging\lan\',[animal 'psych'], 'videoFeatures', [animal num2str(days(day_i)) 'whisk.mat']), ...
                'wsk_time_trace');
            clear 'channels_data';
            clear 'timing';
            
            
            
            blinkfile = fullfile('X:\Lan\Meso-imaging\', animal, [animal '_D' ...
                num2str(days(day_i))  '_blink.mat']);
            if ~exist(blinkfile, 'file')
                disp('no wskile');
                continue;
            end
            blink = load(blinkfile);
            for T=(trialsinds)
                L = find(isnan(blink.blinktrace(:,T)),1);
                x = blink.blinktrace(1:L-1,T);
                blink_time_trace.blinktrace(1,:,T) = interp1(linspace(-2,3,L-1), x, linspace(-2,3,150));
            end
            blink_time_trace.t = linspace(-2, 3, 150);
            
            
            
            save(fullfile('X:\Hadas\Meso-imaging\lan\',[animal 'psych'], 'videoFeatures', [animal num2str(days(day_i)) 'blink_time_trace.mat']), ...
                'blink_time_trace');
        end
        
    
    
end
