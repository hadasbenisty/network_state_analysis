function spike2processing_crispr
addpath(genpath('../spike2_utils'));
addpath('..\functions');
addpath(genpath('../utils'));
addpath(genpath('../meta_data_processing/'));
animals_db = get_animals_meta_data_by_csv;
outputFolder='X:\Hadas\Meso-imaging\CRISPR\traces_data';
[params,channels] = get_channels_param('Vis','F238','blueuv');
CEDS64LoadLib('../spike2_utils/CEDS64ML/');
validsessions = animals_db.isimagingood_list'==find(strcmp(animals_db.isimagingood_lut,'GOOD'));
validsessions=validsessions'&animals_db.toinclude_list==find(strcmp(animals_db.toinclude_lut,'Good'));
for i=1:length(validsessions)
    if validsessions(i) == 0
        continue;
    end
    pre_Session=char(animals_db.folder_list(i));
    if contains(pre_Session,'Control')
        Session = strrep(pre_Session,'9_C','9 C');
        %Session=regexprep(sessions_it{i,:}, '_', ' ')
    elseif contains(pre_Session,'Lhx6')
        Session = strrep(pre_Session,'6_C','6 C');
    elseif contains(pre_Session,'Emx')
        Session = strrep(pre_Session,'x_C','x C');
    end
    tiffsPath=char(fullfile('x:\CardinLab\Antara\',Session));
    disp(strcat('Processing ',Session));

    %% get tiffs, smrx , visual csv path and run the main function
    smrxfilelist = (dir(fullfile(tiffsPath, '*.smrx')));
    dataSmrxFile=fullfile(tiffsPath,smrxfilelist.name);
    outputPath=fullfile(outputFolder,animals_db.folder_list{i});
    mkNewDir(outputPath);
    if ~exist(dataSmrxFile, 'file')||isempty(smrxfilelist)
        disp('no file');
    else
        tic


        process_spike2_crispr('X:\Hadas\Meso-imaging\CRISPR\traces_data\',tiffsPath, dataSmrxFile, outputPath, channels,params,Session);
        toc
        %[timing, channels_data] = spike2_retmappingandwheel(cedpath, outputpath,filename, fsspike2,channels,minRunDuration,minSitDuration,ITITime,pupilSR,tiffsPath,Session);
        %save here
        %spike2_retmapping(filename,Session)
    end
end