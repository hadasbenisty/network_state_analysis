function [tform,R,C,h1,h2] = registerGreenBlue_transform(moving,fixed,R,C)
%register two color images acquired with two cameras 
[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 2000;
metric.NumberOfHistogramBins = 10;
metric.UseAllPixels = false;
optimizer.GrowthFactor = 1.02000;
optimizer.Epsilon = 1.50000e-06;
optimizer.InitialRadius = 6.25000e-03;
optimizer.MaximumIterations = 1000;
se=strel('sphere',6);
minval=min(fixed,[],2);
maxval = max(fixed,[],2);
minval2=min(moving,[],2);
maxval2 = max(moving,[],2);

diffframe1 = reshape(maxval-minval, R, C);
diffframe2 = reshape(maxval2-minval2, R, C);

tform= imregtform(diffframe2,diffframe1,'rigid',optimizer,metric);
toalign_new=imwarp(diffframe2,tform,'OutputView',imref2d(size(diffframe1)),'Fillvalues',0);
frame_new=imwarp(reshape(maxval2, R,C),tform,'OutputView',imref2d(size(diffframe1)),'Fillvalues',0);

h1=figure;subplot(2,2,1),imshow(imadjust(mat2gray(reshape(maxval,R,C))));title('Before_Fixed'); 
subplot(2,2,2),imshow(imadjust(mat2gray(reshape(maxval2,R,C))));title('Before_Moving'); 
subplot(2,2,3),imshow(imadjust(mat2gray(reshape(maxval,R,C))));title('After_Fixed'); 
subplot(2,2,4),imshow(imadjust(mat2gray(frame_new)));title('After_Moving'); 

h2=figure;subplot(2,1,1); imshowpair(imadjust(mat2gray(frame_new)),imadjust(mat2gray(reshape(maxval,R,C))),'Scaling','joint');title('After');
subplot(2,1,2); imshowpair(imadjust(mat2gray(reshape(maxval2,R,C))),imadjust(mat2gray(reshape(maxval,R,C))),'Scaling','joint');title('Before');

end




