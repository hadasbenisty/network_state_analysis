function normalize_pupil_by_dark_light_changes(outputpath, inputpath,overwrite)
addpath(genpath('../pre_processing_scripts/utils'));


output_time_stamp_filename=fullfile(outputpath);
% if contains(output_time_stamp_filename,'Control')
%     Session = strrep(output_time_stamp_filename,'9 C','9_C');
% elseif contains(output_time_stamp_filename,'Lhx6')
%     Session = strrep(pre_Session,'6 C','6_C');
% elseif contains(output_time_stamp_filename,'Emx')
%     Session = strrep(output_time_stamp_filename,'x C','x_C');
% end
display(strcat('pupilsize normalization: ',output_time_stamp_filename));

if exist(fullfile(output_time_stamp_filename,'\roi.mat'), 'file') &&...
        exist(fullfile(output_time_stamp_filename,'\pupil.mat'), 'file') &&...
        ~overwrite
    load(fullfile(output_time_stamp_filename,'\pupil.mat'));
    load(fullfile(output_time_stamp_filename,'\axisroi.mat'));
    %load(strcat(data_time_stamp_filename,'.mat')); % load in the diode information
    %load(strcat(data_time_stamp_filename,'binary.mat'));
    load(fullfile(output_time_stamp_filename,'\roi.mat'))
    t=quantile(roiint,0.95);
    ind=roiint>t | centerx<quantile(centerx,0.01)| centerx>quantile(centerx,0.99)|centery<quantile(centery,0.01)| centery>quantile(centery,0.99);
    areaii(ind)=NaN;
    centerx(ind)=NaN;
    centery(ind)=NaN;
    circularity(ind)=NaN;
    
    areaii_new=medfilt1(areaii,3);
    areaii_new=fillmissing(areaii_new,'nearest');
    areaii_new(areaii_new<quantile(areaii_new,0.01))=NaN;
    areaii_new=fillmissing(areaii_new,'nearest');
    %     doesnt make sense to do since we dont have vis
    %     trace=areaii_new(2000:end);
    %     bl = nanmedian(trace(trace <= quantile(trace, 0.01)));
    %     blmax=nanmedian(trace(trace >= quantile(trace, 0.99)));
    %     areaii_blnorm=(areaii_new-bl)./(blmax-bl);
    %     areaii_bldivide=areaii_new./bl;
    %     areaii_eyenorm=areaii_new./eyesize;
    
    areaii=areaii_new;
    
    save(fullfile(output_time_stamp_filename,'pupil_clean.mat'),'areaii','centerx','centery','circularity');
    disp('saved')
else
    disp('one condition not met')
end
%save(strcat(data_time_stamp_filename,'_normpupil.mat'),'areaii','bl','blmax','areaii_eyenorm','areaii_blnorm','areaii_bldivide');
end

