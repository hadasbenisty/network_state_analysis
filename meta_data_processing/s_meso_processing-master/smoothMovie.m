function[finvideo]= smoothMovie(finvideo,R,C,sigma)
%spatially smooth with a gaussian window
smoothMovie=zeros(R,C,size(finvideo,3));
for t=1:size(finvideo,3)
    currImage=finvideo(:,:,t);
    smoothMovie(:,:,t)=imgaussfilt(currImage,sigma);
end
finvideo=smoothMovie;
end