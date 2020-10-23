function [CondColors,colortouse] = colorpergroup(name)
CondColors=[0.9961,0.5469,0;0,0,0.8008;0.6953,0.1328,0.1328;0.9961,0.8398,0];
if strcmp(name,'highpup')
    colortouse=[0.9961,0.5469,0];
elseif strcmp(name,'lowpup')
    colortouse=[0,0,0.8008];
elseif strcmp(name,'run')
    colortouse=[0.6953,0.1328,0.1328];
elseif strcmp(name,'notrun')
    colortouse=[0.9961,0.8398,0];
end
end

