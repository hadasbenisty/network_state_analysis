function [t, X] = extract_imaging_per_state(stateOn,stateOff,imaging_time, dFoF_parcells)
%% inputs
%stateOn: onset times for state1 (such as locomotion)
%stateOff: offset times for state1
%imaging_time: imaging timestamps 
%dFoF_parcells: imaging data by parcells


X=[];t=[];
for r=1:length(stateOn)
    X=cat(2, X, dFoF_parcells(:,imaging_time>stateOn(r) & imaging_time<stateOff(r)));
    t = cat(1, t, imaging_time(imaging_time>stateOn(r) & imaging_time<stateOff(r)));
end
        
      
      
