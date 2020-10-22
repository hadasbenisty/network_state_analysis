function [Pupil_HighArousal_On_final_loc, Pupil_HighArousal_Off_final_loc,...
    Pupil_LowArousal_On_final_loc, Pupil_LowArousal_Off_final_loc] = ...
    getHighLowonOf4Locomotion(wheelOn_final, wheelOff_final, t_imaging, ...
    Pupil_HighArousal_OnT_int, Pupil_HighArousal_OffT_int, ...
    Pupil_LowArousal_OnT_int, Pupil_LowArousal_OffT_int)


wheelstate = zeros(size(t_imaging));
    for rj=1:length(wheelOn_final)
        wheelstate(findClosestDouble(wheelOn_final(rj), t_imaging):findClosestDouble(wheelOff_final(rj), t_imaging))=1;
    end
    Pupil_HighArousalstate = zeros(size(t_imaging));
    for rj=1:length(Pupil_HighArousal_OnT_int)
        Pupil_HighArousalstate(findClosestDouble(Pupil_HighArousal_OnT_int(rj), t_imaging):findClosestDouble(Pupil_HighArousal_OffT_int(rj), t_imaging))=1;
    end
    high_pup_loc = wheelstate & Pupil_HighArousalstate;
    
    Pupil_LowArousalstate = zeros(size(t_imaging));
    for rj=1:length(Pupil_LowArousal_OnT_int)
        Pupil_LowArousalstate(findClosestDouble(Pupil_LowArousal_OnT_int(rj), t_imaging):findClosestDouble(Pupil_LowArousal_OffT_int(rj), t_imaging))=1;
    end
    low_pup_loc = wheelstate & Pupil_LowArousalstate;
    [Pupil_HighArousal_On_final_loc,Pupil_HighArousal_Off_final_loc]=squaredetect(high_pup_loc,.5);
    [Pupil_LowArousal_On_final_loc,Pupil_LowArousal_Off_final_loc]=squaredetect(low_pup_loc,.5);
    
    Pupil_HighArousal_On_final_loc=t_imaging(Pupil_HighArousal_On_final_loc);
    Pupil_HighArousal_Off_final_loc=t_imaging(Pupil_HighArousal_Off_final_loc);
    
    Pupil_LowArousal_On_final_loc=t_imaging(Pupil_LowArousal_On_final_loc);
    Pupil_LowArousal_Off_final_loc=t_imaging(Pupil_LowArousal_Off_final_loc);
   