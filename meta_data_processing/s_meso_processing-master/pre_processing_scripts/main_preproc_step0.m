%% script for pupil detection and blink amp with nidaq setup
% input file - an excel file specifying which animal on which days should
% be processed. Excel file - Sheet1 should be:
% col A a list of animals and col B a list of days. See 
% \\128.36.220.173\vivo2\Lan\Meso-imaging\lt_tocleanmeso.xlsx for example
% The script analyzes all days specified per animal
clear;
addpath(genpath('utils'));

%% input
overwrite_smrx2binary                                   = false;
overwrite_avimp4_2axisroi                               = false;
overwrite_avimp4__axisroi_2_roi_pupil_whisker           = false;
overwrite_binary2_blink_blinksummary                    = false;
overwirte_binary_roi_2_spontblink                       = false;
overwirte_pupil_axisroi_channels_binary_2_whisker_clean = false;

% smrx channels specifications
ANALOG_RATE = 5000;
CHANNELS_NUM = 7;
channels.BLUE = 1;
channels.UV = 2;
channels.FRAMETICKES = 3;
channels.PHOTO_DIODE = 4;
channels.WHEEL = 5;
channels.AIR_PUFF = 6;
channels.PUPIL_CAMERA = 7;

 
rf=0;
datapath='\\128.36.220.173\vivo2\Lan\Meso-imaging\';% path to imaging + smrx
outputpath = '\\128.36.220.173\vivo2\Hadas\Meso-imaging';% path to save processed mat files
excelfilename = fullfile(datapath,'lt_tocleanmeso.xlsx');% list of animals and days to process


% parameters for extracting wheel motion
MIN_MAX_THRESHOLD = 3;
VOLTS_RANGE = 5;
WINDOW = 2500;
WHEEL_DIAMETER = 15; % cm
WHEEL_DIRECTION = 1;
MAX_WHEEl_SPEED = 60;
SQUARE_DETECT_PARAM = 0.05;
SHORT_BREAK_OR_SHORT_RUN_PARAM = 2500;


%load SON library
cedpath = '../CEDMATLAB/CEDS6ML/';
CEDS64LoadLib(cedpath);
% load animals and days
[daysession,animal]=xlsread(excelfilename);
animal=char(animal);


%% 1. smrx->_binary.mat
% Clean signal matrix
clean_smrx_signal(datapath, outputpath, daysession,animal, channels, rf, ...
            overwrite_smrx2binary, CHANNELS_NUM, WHEEL_DIRECTION, MIN_MAX_THRESHOLD, WINDOW, WHEEL_DIAMETER, MAX_WHEEl_SPEED, ...
            SQUARE_DETECT_PARAM, SHORT_BREAK_OR_SHORT_RUN_PARAM);


%% 2. avi/mp4->_axisroi.mat
% detect axis for pupil
detect_axis_for_pupil(outputpath, datapath, daysession, animal, overwrite_avimp4_2axisroi);

%% 3.  mp4/avi + _axisroi.mat -> _roi.mat, _pupil.mat, _whisker.mat
% extract roi intensity detect pupil size_and_whiskering
extract_roi_intensity_detect_pupil_size_and_whiskering(animal, daysession, ...
    datapath, outputpath, overwrite_avimp4__axisroi_2_roi_pupil_whisker)

%% 4.  _binary.mat + -> _blink + _blinksummary.mat
% eyeblinkdetect-semi automatic
semi_auto_eye_blink(outputpath, daysession, animal, ANALOG_RATE, overwrite_binary2_blink_blinksummary);
%% 5. _binary + _roi -> _spontblink.mat
% spontaneous eyeblink detect, binarize whisking and normalize pupil size
% detects the timing of spont blink
detect_timing_spont_blink(animal, daysession, outputpath, overwirte_binary_roi_2_spontblink);

%% 6. _pupil + _axisroi + channels + _binary ->_whisker_clean
% normalize pupil size based on dark-light changes/minimum pupil/eye size
normalize_pupil_by_dark_light_changes(outputpath, daysession, animal, ...
    overwirte_pupil_axisroi_channels_binary_2_whisker_clean)

