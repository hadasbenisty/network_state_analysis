%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads smr data quired by spike 2. Reads smr files for all animals and
% days as specified in input and saves and file with '_binary' suffix
%
% input - 
%   datapath                       - path to read input files
%   outputpath                     - path to save output files
%   daysession                     - a list of days for processing
%   animal                         - a list of animals as a char array
%   channels                       - a struck specifying channels names and
%                                    numbre. For example 
%                                    channels.BLUE = 1, channels.UV = 2, etc.
%   rf                             - need to ask Lan what this is 
%   overwrite                      - if to overwrite output if exists
%   WHEEL_DIRECTION                - direction of motion of wheel
%   MIN_MAX_THRESHOLD              - need to ask Lan what this is
%   WINDOW                         - need to ask Lan what this is
%   WHEEL_DIAMETER                 - size of wheel
%   MAX_WHEEl_SPEED                - need to ask Lan what this is
%   SQUARE_DETECT_PARAM            - need to ask Lan what this is
%   SHORT_BREAK_OR_SHORT_RUN_PARAM - need to ask Lan what this is
%
%   Written by Lan Tang
%   editted by Hadas Ben Esti 6/28/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function clean_smrx_signal(datapath, outputpath, daysession,animal, channels, rf, ...
            overwrite, CHANNELS_NUM, WHEEL_DIRECTION, MIN_MAX_THRESHOLD, WINDOW, WHEEL_DIAMETER, MAX_WHEEl_SPEED, ...
            SQUARE_DETECT_PARAM, SHORT_BREAK_OR_SHORT_RUN_PARAM)
for filenum=1:size(daysession,1)    
    current_animal=animal(filenum,:);
    data_mat_time_stamp_filename=strcat(current_animal,'_D',num2str(daysession(filenum,1)));
    if ~exist(fullfile(outputpath,current_animal,strcat(data_mat_time_stamp_filename,'_binary.mat')), 'file')...
            || overwrite
        data_smr_time_stamp_filename =fullfile(current_animal, strcat(current_animal,'_D',num2str(daysession(filenum,1))));
        display(strcat('cleaning signal matrix: ',data_smr_time_stamp_filename));
        
            data = process_spike2_smr2mat(datapath, outputpath, data_smr_time_stamp_filename, CHANNELS_NUM);
             
        channels_data = get_channels_data_from_samples(data, channels, rf, ...
            WHEEL_DIRECTION, MIN_MAX_THRESHOLD, WINDOW, WHEEL_DIAMETER, MAX_WHEEl_SPEED, ...
            SQUARE_DETECT_PARAM, SHORT_BREAK_OR_SHORT_RUN_PARAM);
        % Other modalities
        channels_data.startsig = zeros(size(channels_data.diode));
        channels_data.startsig(1:end-500)=channels_data.diode(501:end);        
        channels_data.sound = zeros(size(channels_data.diode));  % no sound
        channels_data.led=zeros(size(channels_data.diode));
        channels_data.CaFrameData =zeros(size(channels_data.diode));
        save(fullfile(outputpath,current_animal,strcat(data_mat_time_stamp_filename,'_binary.mat')),'channels_data');
    end    
end