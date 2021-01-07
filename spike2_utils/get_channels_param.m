function [params,channels] = get_channels_param(sessid,location,color)
                        %% user-selected input parameters
                        params.Loc=location; %recording location, 'B449' or 'F238'
                        params.fsspike2=5000;
                        params.fsimaging=10;
                        params.pupilSR=10;
                        params.batchSize = 3000;
                        params.batches2load=10000;
                        params.regressUV=1; %option of regressing every pixel by UV
                        params.moveLocal=1;%if you want to move the tiff file to a local drive before reading
                        params.drawPlots=0;%plot all event channels and imaging channel, for test purposes
                        params.makeMovies=0;%make tif df/f movies for a segement of the session
                        params.ImageFormat='tif';%choose between cxd or tif
                        % visusal stim
                        params.visStimAn=0; % 1 if vis stim are presented, 0 if not
                        params.airpuffAn=0; % if airpuffs or electrical stim are presented, indicate 1;
                        if strcmp(sessid,'Vis'), params.visStimAn=1;end
                        if strcmp(sessid,'Air'), params.airpuffAn=1;end
                        params.visDur=[];
                        params.visITI=[];
                        %grabs/rcamp vs grabs only
                        params.signalsExtraction.firstCaFrame = 1;
                        params.signalsExtraction.sigs = color;% 'blueuv' or 'RCaMP_AC'
                        params.signalsExtraction.blueuvRatio = 1;
                        params.resizeRatio = 0.5;
                        % detrending and spatial filtering params
                        params.deterend.filtLen = 150;
                        params.deterend.filtcutoff = 0.001;
                        params.deterend.method = 'FIR';
                        params.spatiotemporalfilter=1; %set it to 1 if you want to do spatial filtering  otherwise 0
                        
                        %run params
                        params.minRunDuration=5; %in seconds
                        params.minSitDuration=10; %in seconds
                        params.ITITime=5; %  dead time after events in seconds (so they dont overlap)
                        
                        %analysis window for events (airpuff and visstim)
                        params.preEventWin=5;
                        params.postEventWin=5;
                        params.baselineWin=2; %use this period to calculate baseline during within trial normalization
                        params.angle2rotate=90;
                        % spike2 channels for cardin 238
                        channels.BLUE = 1;%blue led
                        channels.UV = 2;%uv led
                        channels.FRAMETICKS = 3;%green/uv mesocam triggers
                        channels.PHOTO_DIODE = 4;%visual stim trigger
                        channels.WHEEL = 5;
                        channels.AIR_PUFF = 6; %this is either airpuff channel or electrical stim channel
                        channels.PUPIL_CAMERA = 7;
end

