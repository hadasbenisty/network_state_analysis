function [tform] = registerGreenBlue(moving,fixed,R,C,outputPath)
%register two color images acquired with two cameras 
minval=min(fixed,[],2);
maxval = max(fixed,[],2);
minval2=min(moving,[],2);
maxval2 = max(moving,[],2);

fixed_diffframe1 = reshape(maxval-minval, R, C);
moving_diffframe2 = reshape(maxval2-minval2, R, C);

[MOVINGREG] = registerImages(moving_diffframe2,fixed_diffframe1);

tform=MOVINGREG.Transformation; 
moving_reg=MOVINGREG.RegisteredImage;

%create a flickering tif file to see before and after registration 
h1=mat2gray(fixed_diffframe1); 
outFile=fullfile(outputPath,'Before_reg.tif');
imwrite(im2uint16(h1),outFile ,'tif');
h1 = mat2gray(moving_diffframe2); 
imwrite(im2uint16(h1), outFile,'WriteMode', 'append');

h2=mat2gray(fixed_diffframe1); 
outFile=fullfile(outputPath,'After_reg.tif');
imwrite(im2uint16(h2),outFile ,'tif');
h2 = mat2gray(moving_reg); 
imwrite(im2uint16(h2), outFile,'WriteMode', 'append');

h3=figure;subplot(2,2,1),imshow(imadjust(mat2gray(reshape(maxval,R,C))));title('Before Fixed'); 
subplot(2,2,2),imshow(imadjust(mat2gray(reshape(maxval2,R,C))));title('Before Moving'); 
subplot(2,2,3),imshow(imadjust(mat2gray(reshape(maxval,R,C))));title('After Fixed'); 
subplot(2,2,4),imshow(imadjust(mat2gray(moving_reg)));title('After Moving'); 
outFile=fullfile(outputPath,'GreenBlueRegistration.fig');
saveas(h3,outFile,'fig'); 
end




