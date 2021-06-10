function [tform, blue_mov_trans, uv_mov_trans] = parcellateData(blue_mov, uv_mov)
[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 2000;
metric.NumberOfHistogramBins = 10;
metric.UseAllPixels = false;
optimizer.GrowthFactor = 1.02000;
optimizer.Epsilon = 1.50000e-06;
optimizer.InitialRadius = 6.25000e-03;
optimizer.MaximumIterations = 1000;
se=strel('sphere',6);
minblue=min(blue_mov,[],3);
maxblue = max(blue_mov,[],3);
% first use whole frame to registrate
toalign= double(imadjust(imdilate(imerode(maxblue-minblue,se),se)));
% toalign2=double(im2bw(toalign,0.1));
template=imresize(imread('cortex_template.tif'),0.5,'bilinear');
template2=double(im2bw(template,0.1));
figsize=size(template);
tform= imregtform(toalign,template2,'similarity',optimizer,metric);

toalign_new=imwarp(toalign,tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
blue_new=imwarp(maxblue,tform,'OutputView',imref2d(size(template)),'Fillvalues',0);

h=figure('Position',[50,50,1200,800]);
subplot(2,3,1);
imshow(imadjust(toalign_new));
title('diffblue');
subplot(2,3,4);
imshowpair(imadjust(toalign_new),uint8(5*template2));
subplot(2,3,2);
imshow(imadjust(blue_new));
title('blue');
subplot(2,3,5);
imshowpair(imadjust(blue_new),uint8(5*template2));
if exist('uv_mov','var')&&~isempty(uv_mov)
    firstuvframe = uv_mov(:,:,1);
    uv_new=imwarp(firstuvframe,tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
    
    subplot(2,3,3);
    imshow(imadjust(uv_new));
    title('uv_mov');
    subplot(2,3,6);
    imshowpair(imadjust(uv_new),uint8(5*template2));
end
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
        imshow(imadjust(maxblue-minblue),'InitialMagnification','fit');
        drawmask= drawpolygon(gca);
        %mask=find(reshape(drawmask.createMask,[],1));
        mask=drawmask.createMask;
        close(h2);
        aframe=maxblue;
%         bframe=minblue;
        aframe(mask==0)=0;
%         bframe(mask==0)=0;
        toalign= imadjust(imdilate(imerode(aframe,se),se));
%         toalign2=double(im2bw(toalign,0.1));
        template3=template2;
        if nanmedian(find(nanmean(single(mask))))>figsize(1)/3*2
            template3(:,1:round(figsize(1)/2))=0;
        elseif nanmedian(find(nanmean(single(mask))))<figsize(1)/3
            template3(:,round(figsize(1)/2):end)=0;
        end
        
        tf_T=0;
        for ireg=1:5
            tform= imregtform(toalign,template3,'similarity',optimizer,metric);
            tf_T=tf_T+tform.T/5;
        end
        tform.T=round(tf_T,2);
        
        toalign_new=imwarp(imadjust(imdilate(imerode(maxblue-minblue,se),se)),tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
        blue_new=imwarp(maxblue,tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
        
        figure('Position',[50,50,1200,800]);
        subplot(2,3,1);
        imshow(imadjust(toalign_new));
        title('diffblue');
        subplot(2,3,4);
        imshowpair(imadjust(toalign_new),uint8(5*template2));
        subplot(2,3,2);
        imshow(imadjust(blue_new));
        title('blue');
        subplot(2,3,5);
        imshowpair(imadjust(blue_new),uint8(5*template2));
        if exist('uv_mov','var')&&~isempty(uv_mov)
            
            uv_new=imwarp(firstuvframe,tform,'OutputView',imref2d(size(template)),'Fillvalues',0);
            
            subplot(2,3,3);
            imshow(imadjust(uv_new));
            title('uv_mov');
            subplot(2,3,6);
            imshowpair(imadjust(uv_new),uint8(5*template2));
        end
        
end
if exist('uv_mov','var')&&~isempty(uv_mov)
    
    
    uv_mov_trans = transform_frames(uv_mov, tform);
else
    uv_mov_trans=[];
end
blue_mov_trans = transform_frames(blue_mov, tform);

end



