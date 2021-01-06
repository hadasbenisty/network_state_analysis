function main_time_traces_to_spont_2p_lan


addpath(genpath('../meta_data_processing/'));

T = readtable('X:\Lan\FMB208 conditioning imaging data\animal_log_imaging_lt.xlsx');
animalNames = T.AnimalID;

fsspike2=5e3;


for animal_i =1:length(animalNames)
    
    fsimaing=10;
    delay_filt=0;
    
    animal = animalNames{animal_i};
    dayslist = eval(T.PsychTest{animal_i});
    for day_i = 1:length(dayslist)
        disp([animalNames{animal_i}  ' '  num2str(dayslist(day_i))]);
        
        datapath = ['X:\Hadas\Meso-imaging\lan\' animalNames{animal_i} 'psych\spt' ];
        
        if exist(fullfile(datapath, [animalNames{animal_i} num2str(dayslist(day_i)) 'imaging_time_traces_global_ITI_dfff.mat']), 'file')
            d=load(fullfile(datapath, [animalNames{animal_i} num2str(dayslist(day_i)) 'imaging_time_traces_global_ITI_dfff.mat']));
            if isfield(d.imaging_time_traces,'grid4')
                continue;
            end
        end
        imagingfile = fullfile(datapath, [animalNames{animal_i} '_D' num2str(dayslist(day_i)) '_full_dfff.mat']);
       
        
        if ~exist(imagingfile, 'file')
            disp('files missing');
            continue;
        end
        
        
        % spike2
        spike2pth = fullfile('X:\Hadas\Meso-imaging\lan\spike2data', animalNames{animal_i});

        if ~exist(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(dayslist(day_i)) '.mat']),'file')
            
            
                disp('No smrx');
                continue;
           
        else
            load(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(dayslist(day_i)) '.mat']),'channels_data',...
                'timing', 't_imaging');
            
        end
        before_win = 1;
        after_win = 3;
        
        data = load(imagingfile);
        
        if length(t_imaging) > size(data.imagingtraces,2)
            t_imaging=t_imaging(1:size(data.imagingtraces,2));
        end
        X =data.imagingtraces(:, 1:length(t_imaging));
       
       
        stim_timestamps = timing.stimstart/fsspike2;
        stim_timestamps=stim_timestamps(1:75);
        [t, imaging_time_traces.cells]                       = time_trace2ITI(X, t_imaging, stim_timestamps, before_win, after_win, fsimaing);
       
        
      
        save(fullfile(datapath, [animalNames{animal_i} num2str(dayslist(day_i)) 'imaging_time_traces_global_ITI_dfff.mat']),...
            'imaging_time_traces' , 't');
        clear imaging_time_traces;
        
        continue;
        t_spike2 = [1:length(channels_data.wheel)]';
        [ond,ofd] = squaredetect(channels_data.diode,.5);
        
        [~, spike2time_traces.wheelspeed] = time_trace2ITI(channels_data.wheelspeed.', t_spike2, ond, before_win, after_win, fsspike2);
        [spike2time_traces.t, spike2time_traces.wheel] = time_trace2ITI(channels_data.wheel.', t_spike2, ond, before_win, after_win, fsspike2);
        spike2time_traces.t=spike2time_traces.t/fsspike2;
        mkdir(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features'));
        save(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features', [animalNames{animal_i} num2str(dayslist(day_i)) 'spike2_ITI.mat']),  'spike2time_traces');
        clear spike2time_traces;
        
        
        % whisking
        wskile = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(dayslist(day_i))  '_whisker_clean.mat']);
        if ~exist(wskile, 'file')
            disp('no wskile');
            continue;
        end
        wsk = load(wskile,  'whisk_binary');
        wsk.whisk_binary =wsk.whisk_binary(1:min(length(wsk.whisk_binary),length(t_imaging)));
        [~, wsk_time_trace.binary] = time_trace2ITI(wsk.whisk_binary', t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        
        wskile = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(dayslist(day_i))  '_whisker.mat']);
        if ~exist(wskile, 'file')
            disp('no wskile');
            continue;
        end
        wsk = load(wskile,  'whisker');
        wsk.whisker =wsk.whisker(1:min(length(wsk.whisker),length(t_imaging)));
        
        [wsk_time_trace.t, wsk_time_trace.whisker] = time_trace2ITI(wsk.whisker', t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        
        mkdir(fullfile('X:\Hadas\Meso-imaging\lan\',[animalNames{animal_i} 'psych'], 'videoFeatures'))
        save(fullfile('X:\Hadas\Meso-imaging\lan\',[animalNames{animal_i} 'psych'], 'videoFeatures', [animalNames{animal_i} num2str(dayslist(day_i)) 'whisk_ITI.mat']), ...
            'wsk_time_trace');
        clear 'channels_data';
        clear 'timing';
        
        
        
        
    end
    
    
    
end
