function [spike2points,caframes,pupframes]=getandalign_lt(t1,t2,spike2timevec,catimevec,puptvec)
if ~isempty(catimevec)
spike2points=find(spike2timevec>=t1 & spike2timevec<=t2); 
caframes=find(catimevec>=t1 & catimevec<=t2);
pupframes=find(puptvec>=t1 & puptvec<=t2);
else
spike2points=find(spike2timevec>=t1 & spike2timevec<=t2); 
caframes=[];
pupframes=find(puptvec>=t1 & puptvec<=t2); 
end





