function [tform,R,C] = get_alignment_transform(insig,params)
[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 2000;
metric.NumberOfHistogramBins = 10;
metric.UseAllPixels = false;
optimizer.GrowthFactor = 1.02000;
optimizer.Epsilon = 1.50000e-06;
optimizer.InitialRadius = 6.25000e-03;
optimizer.MaximumIterations = 1000;
se=strel('sphere',6);
minval=min(insig,[],2);
maxval = max(insig,[],2);
template=imresize(imread('cortex_template.tif'),0.5,'bilinear');

[R,C] = size(template);
if strcmp(params.dff,'baselineFilt')
diffframe = reshape(maxval-minval, R, C);
elseif strcmp(params.dff,'min10%')
diffframe = reshape(mean(insig,2), R, C);
end 
% first use whole frame to registrate
toalign= double(imadjust(imdilate(imerode(diffframe,se),se)));
% toalign2=double(im2bw(toalign,0.1));
template2=double(im2bw(template,0.1));
tform= imregtform(toalign,template2,'similarity',optimizer,metric);

toalign_new=imwarp(toalign,tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
frame_new=imwarp(reshape(maxval, R,C),tform,'OutputView',imref2d(size(template)),'Fillvalues',0);

h=figure('Position',[50,50,1200,800]);
subplot(2,2,1);
imshow(imadjust(toalign_new));
title('diff');
subplot(2,2,3);
imshowpair(imadjust(toalign_new),uint8(5*template2));
subplot(2,2,2);
imshow(imadjust(frame_new));
title('sig');
subplot(2,2,4);
imshowpair(imadjust(frame_new),uint8(5*template2));

ButtonName = questdlg('Is automatic OK?', ...
    '', ...
    'Yes','No','Cancel','Yes');
% Handle response
switch ButtonName
    case 'Cancel'
        return;
    case 'Yes'
        % do nothing;
    case 'No'
        close(h);
        % create a circular mask, anything outside of the mask is set to black
        h2=figure('Position',[300,300,700,700]);
        imshow(imadjust(diffframe),'InitialMagnification','fit');
        drawmask= drawpolygon(gca);
        %mask=find(reshape(drawmask.createMask,[],1));
        mask=drawmask.createMask;
        close(h2);
        aframe=reshape(maxval, R, C);
%         bframe=minval;
        aframe(mask==0)=0;
%         bframe(mask==0)=0;
        toalign= imadjust(imdilate(imerode(aframe,se),se));
%         toalign2=double(im2bw(toalign,0.1));
        template3=template2;
        if nanmedian(find(nanmean(single(mask))))>R/3*2
            template3(:,1:round(R/2))=0;
        elseif nanmedian(find(nanmean(single(mask))))<R/3
            template3(:,round(R/2):end)=0;
        end
        
        tf_T=0;
        for ireg=1:5
            tform= imregtform(toalign,template3,'similarity',optimizer,metric);
            tf_T=tf_T+tform.T/5;
        end
        tform.T=round(tf_T,2);
        
        toalign_new=imwarp(imadjust(imdilate(imerode(diffframe,se),se)),tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
        frame_new=imwarp(reshape(maxval, R,C),tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
        
        figure('Position',[50,50,1200,800]);
        subplot(2,2,1);
        imshow(imadjust(toalign_new));
        title('diffblue');
        subplot(2,2,4);
        imshowpair(imadjust(toalign_new),uint8(5*template2));
        subplot(2,2,2);
        imshow(imadjust(frame_new));
        title('sig');
        subplot(2,2,3);
        imshowpair(imadjust(frame_new),uint8(5*template2));
        
        
end




