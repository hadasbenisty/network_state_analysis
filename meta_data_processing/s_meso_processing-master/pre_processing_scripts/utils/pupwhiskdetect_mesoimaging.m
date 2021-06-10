function [roiint,areaii,centerx,centery,circlemetric,whisker,nframes]=pupwhiskdetect_mesoimaging(filename,eyerect,facerect,roi_eye,roi_face,thres,thres2,thres3)

display(strcat('pupil detection :  ',filename));
mov=VideoReader(filename);

pointdist=@(a,b) sqrt(sum((a-b).^2));
nframes=mov.NumberofFrames;
roiint=zeros(nframes,1);
centers=nan(nframes,2);
areaii=nan(nframes,1);
circlemetric=zeros(nframes,1);
whisker=nan(nframes,1);

linecenter=[-eyerect(3)+eyerect(4),-eyerect(1)+eyerect(2)]./2;
edges=-0.5:2:255.5;
warning off;
for iframe=1:nframes
    final=0;
    img=rgb2gray(read(mov,iframe));
    img_eye=img(eyerect(1):eyerect(2),eyerect(3):eyerect(4));
    imgeye=255-imfill(255-img_eye);
    binarizedimg=im2bw(imgeye,thres);
    roiint(iframe)=mean(binarizedimg(roi_eye(eyerect(1):eyerect(2),eyerect(3):eyerect(4))));%measure intensity within ROI
    
    imgpup=1-im2bw(imgeye,thres2);
    [B,L,n] = bwboundaries(imgpup);
    
    stats = regionprops(L,'Area','Centroid'); %Estimate each object's area and perimeter
    index=zeros(n,4); %area, circularity, amd centroid
    % Three rules to validate a pupil object
    % 1. A simple metric indicating the roundness of an object: circularity=4*pi*area/perimeter^2.
    % 2. Area change does not exceed limitation (square-root area ratio <1.1)
    % 3. Centroid change does not exceed limitation (cetroid distance from last fram < 5 pix)
    lastarea=areaii(max(iframe-1,1));
    lastcentroid=centers(max(iframe-1,1),:);
    circularthres=0.4;
    
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
        
        if isnan(lastarea) & area>30
            index(k,1)=area;
        elseif (area/lastarea<((sqrt(area/pi)+10)/sqrt(area/pi))^2) & (area/lastarea>((sqrt(area/pi)-10)/sqrt(area/pi))^2)
            index(k,1)=area;
        end
        
        if circularity > circularthres
            index(k,2)=circularity;
        end
        
        if isnan(lastcentroid(1)) & centroid(1)<linecenter(1)+50 & centroid(2)>linecenter(2)-50 &centroid(2)<linecenter(2)+50 & pointdist(linecenter,centroid)<50
            index(k,3:4)=centroid;
        elseif    pointdist(lastcentroid,centroid)<15
            index(k,3:4)=centroid;
        end
        
        if sum(index(k,:)>0)==4 & circularity>circlemetric(iframe)
            areaii(iframe)=area;
            centers(iframe,:)=centroid;
            circlemetric(iframe)=circularity;
            final=k;
        end
    end
    
    %          imshow(imgeye);
    %          hold on;
    %          if final>0
    %              bound=B{final};
    %              plot(bound(:,2), bound(:,1), 'r', 'LineWidth', 2);
    %          end
    %          waitforbuttonpress;
    
    
    % calculate frame-to-frame correlation of nose roi
    imgface=img(facerect(1):facerect(2),facerect(3):facerect(4));
    %imgface=imgaussfilt(imgface,4);
    %figure(1);imshow(imadjust(imgface));
    currentface=reshape(im2bw(imgface,thres3),[],1);
    whisker(iframe)=sum(currentface);

    if mod(iframe,1000)==0
        display(strcat('pupil detection ', filename,'  ',num2str(iframe), ' out of ', num2str(nframes)));
    end
    
end

centery=centers(:,2)-linecenter(2);
centerx=centers(:,1)-linecenter(1);
circlemetric(circlemetric==0)=NaN;

clear mov;
