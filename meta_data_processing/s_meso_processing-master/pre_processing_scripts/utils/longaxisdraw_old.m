function [longaxis,roi,imagthres,pupilthres]=longaxisdraw_old(moviename) %.mp4 or .avi
mov=VideoReader(moviename);

%for k=linspace(nframes/2-nframes/4,nframes/2+nframes/4,100)
for k=linspace(100,1100,10)
    img = rgb2gray(read(mov,round(k)));
    %if k<=nframes/2-nframes/4+1
    if k==100
        imag=round(img/10);
    else
        imag=imag+round(img/10);
    end
end
imshow(imag);

button=0; 
while button~=49
    [~,~,button]=ginput(1);
    stframe=1000*rand();
    for k=linspace(stframe,min(stframe+1000,2000),10)
        img = rgb2gray(read(mov,round(k)));
    %if k<=nframes/2-nframes/4+1
        if k==stframe
           imag=round(img/10);
        else
            imag=imag+round(img/10);
        end
    end
imshow(imag);
end


% get the long and short axis of the eye
longaxis=imline();
h=wait(longaxis); %h is the endpoint position of line

%1-2 col: long axis, 3 col: image size,
longaxis=[h,[size(imag,1);size(imag,2)]];

%select ROI manually
h3= impoly(gca);
wait(h3);
roi=h3.createMask;

% set binary threshold
% delete light spot
inroi=imag(roi);
outroi=imag(~roi);

imagthres=mean(inroi)./mean(outroi);

button=0;
if imagthres>=1||imagthres<=0
imagthres=0.7;
else
end
imagbin=im2bw(imag,imagthres);

figure(2);
imshow(imagbin);

while button~=49
    [~,~,button]=ginput(1);
    if button==43
        imagthres=min(imagthres+0.025,1);
        imagbin=im2bw(imag,imagthres);
    elseif button==45
        imagthres=imagthres-0.025;
        imagbin=im2bw(imag,imagthres);
    end
    figure(2);
    imshow(imagbin);
end

imagbin=im2bw(imag,imagthres);
figure(2);
imshow(imagbin);


%now adjust the threshold for pupil detection
figure(3);
button=0;
pupilthres=imagthres/2;
imag = rgb2gray(read(mov,2000));
imag=255-imfill(255-imag);
imag = imgaussfilt(imag,4);

pupilbin=im2bw(imag,pupilthres);
imshow(pupilbin);
while button~=49
    [~,~,button]=ginput(1);
    if button==43
        pupilthres=min(pupilthres+0.005,1);
        pupilbin=im2bw(imag,pupilthres);
    elseif button==45
        pupilthres=pupilthres-0.005;
        pupilbin=im2bw(imag,pupilthres);
    elseif button==50
        imag=rgb2gray(read(mov,round(2000+1000*rand())));
        imag=255-imfill(255-imag);
        imag = imgaussfilt(imag,4);
        pupilbin=im2bw(imag,pupilthres);
    end
    figure(3);
    imshow(pupilbin);
end

pupilbin=im2bw(imag,pupilthres);
figure(3);
imshow(pupilbin);

close all;