function main_time_traces_to_spont


addpath(genpath('../meta_data_processing/'));
% load
[~, allen_parcels] = getParcellsByLansAllansAtlas;


% fltstr = 'pixelwise';
fltstr = 'spt';
%
animalNames = {'xx' 'xt' 'xu'   'xs' 'xw'  'xz' };% };%   };%    };
fsspike2=5e3;

spike2pth0 = 'X:\Hadas\Meso-imaging\lan\spike2data';
% addpath(genpath('../../../neuronal-behavioral-processing/utils/'));
dbstop if error;


data_smr_path = 'X:\Lan\Meso-imaging\';
% cedpath = '../pre_processing_scripts/utils/CEDS64ML';

for animal_i =1:length(animalNames)
    switch animalNames{animal_i}
        case {'xu','xv','xt','xs'}
            fsimaing=33;
            delay_filt = 500;
        otherwise
            fsimaing=10;
            delay_filt=150;
    end
    animal = animalNames{animal_i};
    [~,days] = animaltodays(animal);
    for day_i = 1:length(days)
        disp([animalNames{animal_i}  ' '  num2str(days(day_i))]);
        
        datapath = ['X:\Hadas\Meso-imaging\lan\' animalNames{animal_i} 'psych\' fltstr '\'];
        spike2pth = fullfile(spike2pth0, animalNames{animal_i});
        
        if exist(fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'imaging_time_traces_global_ITI_dfff.mat']), 'file')
            d=load(fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'imaging_time_traces_global_ITI_dfff.mat']));
            if isfield(d.imaging_time_traces,'grid4')
                continue;
            end
        end
        parfile = fullfile(datapath, [animalNames{animal_i} '_' num2str(days(day_i)) '_allen_dfff.mat']);
        parfileGal = fullfile(datapath,'gal',[animalNames{animal_i} '_' num2str(days(day_i)) ...
            '_global_dfff.mat']);
        parfileGrid = fullfile(datapath,[animalNames{animal_i} '_' num2str(days(day_i)) ...
            '_grid4_dfff.mat']);
        
        if ~exist(parfile, 'file')||~exist(parfileGal, 'file')|| ~exist(parfileGrid, 'file')
            disp('files missing');
            continue;
        end
        
        
        % spike2
        if ~exist(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(days(day_i)) '.mat']),'file')
            
            filename = fullfile(data_smr_path, animalNames{animal_i}, [animalNames{animal_i} '_D' num2str(days(day_i))]);
            if ~exist([filename '.smrx'],'file')
                disp('No smrx');
                continue;
            end
            [timing, channels_data] = process_spike2(cedpath, '',filename, fsspike2);
            t_spike2 = linspace(0, length(channels_data.startsig)-1,length(channels_data.startsig))/fsspike2;
            t_imaging = timing.mesostart/fsspike2;
            
            save(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(days(day_i)) '.mat']), 't_spike2', 'channels_data', 't_imaging',...
                'timing', 't_spike2');
        else
            load(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(days(day_i)) '.mat']),'channels_data',...
                'timing', 't_imaging');
            
        end
        before_win = 1;
        after_win = 3;
        
        pardataAllan = load(parfile);
        pardata_gal = load(parfileGal);
        pardata_grid = load(parfileGrid);
        
        [parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
        regionLabel.Allen = allen_parcels.regionNum;
        regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
        % [roiLabelsbyAllen_gal, regionLabel.Gal, maskByAllen_gal, maskByAllen.Gal] = get_gal_parcels_lables(animal);
        [parcels_names.Gal, finalindex.Gal, maskByAllen.Gal, regionLabel.Gal, roiLabelsbyAllen.Gal] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
            regionLabel.Allen, animal);
        
        
        
        
        if fsimaing < 33
            t_imaging=t_imaging(1:2:end);
        end
        t_imaging = t_imaging(1:end-delay_filt);
        
        if length(t_imaging) > size(pardataAllan.parcels_time_trace,2)
            t_imaging=t_imaging(1:size(pardataAllan.parcels_time_trace,2));
        end
        X.Allen =pardataAllan.parcels_time_trace(:, 1:length(t_imaging));
        X.Gal =pardata_gal.parcels_time_trace(:, 1:length(t_imaging));
        X.grid4 =pardata_grid.parcels_time_trace(:, 1:length(t_imaging));
        X.Allen = X.Allen(finalindex.Allen,:);
        X.Gal = X.Gal(finalindex.Gal,:);
        X.grid4 = X.grid4(finalindex.grid4,:);
        
        stim_timestamps = timing.stimstart/fsspike2;
        stim_timestamps=stim_timestamps(1:75);
        [~, imaging_time_traces.grid4]                       = time_trace2ITI(X.grid4, t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        [~, imaging_time_traces.Gal]                       = time_trace2ITI(X.Gal, t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        [imaging_time_traces.t, imaging_time_traces.Allen] = time_trace2ITI(X.Allen, t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        
        
        
        %         parcelsLabels.Allen = allen_parcels.regionNum;
        roiLabelsbyAllen.Allen = 1:length(allen_parcels.names);
        regionLabel.nameslegend = {'Rest','Visual','Parietal','Temp','Aud','R-S','S-S','Motor'};
        
        save(fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'imaging_time_traces_global_ITI_dfff.mat']), 'imaging_time_traces' ,...
            'roiLabelsbyAllen', 'maskByAllen', 'regionLabel');
        clear imaging_time_traces;
        clear parcelsLabels;
        clear trialslabels;
        continue;
        t_spike2 = [1:length(channels_data.wheel)]';
        [ond,ofd] = squaredetect(channels_data.diode,.5);
        
        [~, spike2time_traces.wheelspeed] = time_trace2ITI(channels_data.wheelspeed.', t_spike2, ond, before_win, after_win, fsspike2);
        [spike2time_traces.t, spike2time_traces.wheel] = time_trace2ITI(channels_data.wheel.', t_spike2, ond, before_win, after_win, fsspike2);
        spike2time_traces.t=spike2time_traces.t/fsspike2;
        mkdir(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features'));
        save(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features', [animalNames{animal_i} num2str(days(day_i)) 'spike2_ITI.mat']),  'spike2time_traces');
        clear spike2time_traces;
        
        
        % whisking
        wskile = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(days(day_i))  '_whisker_clean.mat']);
        if ~exist(wskile, 'file')
            disp('no wskile');
            continue;
        end
        wsk = load(wskile,  'whisk_binary');
        wsk.whisk_binary =wsk.whisk_binary(1:min(length(wsk.whisk_binary),length(t_imaging)));
        [~, wsk_time_trace.binary] = time_trace2ITI(wsk.whisk_binary', t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        
        wskile = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(days(day_i))  '_whisker.mat']);
        if ~exist(wskile, 'file')
            disp('no wskile');
            continue;
        end
        wsk = load(wskile,  'whisker');
        wsk.whisker =wsk.whisker(1:min(length(wsk.whisker),length(t_imaging)));
        
        [wsk_time_trace.t, wsk_time_trace.whisker] = time_trace2ITI(wsk.whisker', t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        
        mkdir(fullfile('X:\Hadas\Meso-imaging\lan\',[animalNames{animal_i} 'psych'], 'videoFeatures'))
        save(fullfile('X:\Hadas\Meso-imaging\lan\',[animalNames{animal_i} 'psych'], 'videoFeatures', [animalNames{animal_i} num2str(days(day_i)) 'whisk_ITI.mat']), ...
            'wsk_time_trace');
        clear 'channels_data';
        clear 'timing';
        
        
        
        
    end
    
    
    
end
