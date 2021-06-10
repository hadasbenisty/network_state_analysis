function [roiint,areaii,centerx,centery,circlemetric,framenum]=pupildetect_bw(filename,longaxis,roi,thres,thres2)
display(strcat('pupil detection :  ',filename));
mov=VideoReader(filename);

pointdist=@(a,b) sqrt(sum((a-b).^2));
framenum=mov.NumberofFrames;
roiint=zeros(framenum,1);
centers=nan(framenum,2);
areaii=nan(framenum,1);
circlemetric=zeros(framenum,1);

h=longaxis(:,1:2); %first 2 columns: coordinates of long axis
sizeimg=longaxis(:,3); %third column: image size
linecenter=round((h(1,:)+h(2,:))/2); %long line center

warning off;
for i=1:framenum
    final=0;
    img=rgb2gray(read(mov,i));
    
    imgeye=255-imfill(255-img);
    binarizedimg=im2bw(imgeye,thres);
    roiint(i)=mean(binarizedimg(roi));%measure intensity within ROI
           
    img2 = imgaussfilt(imgeye,4); %denoise
    %img2=imadjust(img2);       
   
    %imgnew2=bwareaopen(img2,30);
    imgnew2=1-im2bw(img2,thres2);
    %se = strel('disk',5);
    %imgnew2 = imclose(imgnew2,se);
    %imgnew2=imdilate(imgnew2,se);
    %imgnew2=imfill(imgnew2,'holes');
    
    %[B,L,n] = bwboundaries(1-img2,'noholes'); % boundary and area info in B and L
    [B,L,n] = bwboundaries(imgnew2);
    
    stats = regionprops(L,'Area','Centroid'); %Estimate each object's area and perimeter 
    index=zeros(n,4); %area, circularity, amd centroid
    % Three rules to validate a pupil object
    % 1. A simple metric indicating the roundness of an object: circularity=4*pi*area/perimeter^2.
    % 2. Area change does not exceed limitation (square-root area ratio <1.1)
    % 3. Centroid change does not exceed limitation (cetroid distance from last fram < 5 pix)
    lastarea=areaii(max(i-1,1));
    lastcentroid=centers(max(i-1,1),:);
    circularthres=0.4; %%%%%%%%%%%%%%% set cirularity threshold

    % loop over the boundaries and determine if potential objects' area, centroid and
    % circularity meet criteria
    for k = 1:n
        % calculate circularity
        boundary = B{k}; % obtain (X,Y) boundary coordinates corresponding to 'k'th object
        delta_sq = diff(boundary).^2;
        perimeter = sum(sqrt(sum(delta_sq,2)));               
        area = stats(k).Area;% obtain the area calculation corresponding to label 'k'        
        centroid = stats(k).Centroid;
        circularity = 4*pi*area/perimeter^2; %roundness metric
        
        if isnan(lastarea) & area>200%%%%%%%%%%%%%% %estimate the size of pupil
            index(k,1)=area;
        elseif (area/lastarea<((sqrt(area/pi)+10)/sqrt(area/pi))^2) & (area/lastarea>((sqrt(area/pi)-10)/sqrt(area/pi))^2)
           index(k,1)=area;
        end
        
        if circularity > circularthres
            index(k,2)=circularity;
        end
        
        if isnan(lastcentroid) & centroid(1)<linecenter(1)+50 & ...
                centroid(2)>linecenter(2)-50 &centroid(2)<linecenter(2)+50 &...
                pointdist(linecenter,centroid)<50 %%%%%%%%%%% pupil within central region
            index(k,3:4)=centroid;
        elseif    pointdist(lastcentroid,centroid)<15 %%%%%%%%%%estimate frame2frame pupil movement
            index(k,3:4)=centroid;  
        end
        
        if sum(index(k,:)>0)==4 & circularity>circlemetric(i)
            areaii(i)=area;
            centers(i,:)=centroid;
            circlemetric(i)=circularity;
            final=k;
        end
       
    end
 
% 
     if mod(i,400)==0
         display(strcat('pupil detection ', filename,'  ',num2str(i), ' out of ', num2str(framenum)));
%         imshow(img2);
%         hold on;
%         if final>0
%             bound=B{final};
%             plot(bound(:,2), bound(:,1), 'r', 'LineWidth', 2);
%         end
%     waitforbuttonpress;    
     end
end

centery=centers(:,2)-linecenter(2);
centerx=centers(:,1)-linecenter(1);

circlemetric(circlemetric==0)=NaN;

clear mov;
