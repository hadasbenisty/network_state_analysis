function [roi_eye,roi_face,eyerect,facerect,eyethres,pupilthres,facethres]=longaxisdraw_meso_new(moviename) %.mp4 or .avi
mov=VideoReader(moviename);
disp('processing')
%for k=linspace(nframes/2-nframes/4,nframes/2+nframes/4,100)
imag = rgb2gray(read(mov,100));
figure(1);
imshow(imag);
button=0;
[~,~,button]=ginput(1);
while button~=49
    stframe=10000*rand();
    imag = rgb2gray(read(mov,round(stframe)));
    f1=figure(1);
    imshow(imag);
    [~,~,button]=ginput(1);
end

% circle fa1ce region
h1= imrect(gca,[770,120,200,100]);
wait(h1);
roi_face=h1.createMask;
cood=round(h1.getPosition);
facerect=[cood(2),cood(2)+cood(4),cood(1),cood(1)+cood(3)];

figure(2);
button=0;
facethres=0.4;
imag_face_s=imag(cood(:,2):cood(:,2)+cood(:,4),cood(:,1):cood(:,1)+cood(:,3));
%imag_face_s=225-imfill(225-imag_face_s,'holes');
facebin=im2bw(imag_face_s,facethres);
subplot(2,1,1);
imshow(imag_face_s);
subplot(2,1,2);
imshow(facebin);
[~,~,button]=ginput(1);

while button~=49
    if button==43
        facethres=min(facethres+0.01,1);
        facebin=im2bw(imag_face_s,facethres);
    elseif button==45
        facethres=facethres-0.01;
        facebin=im2bw(imag_face_s,facethres);
    elseif button==50
        imag=rgb2gray(read(mov,round(2000+1000*rand())));
        imag_face_s=imag(cood(:,2):cood(:,2)+cood(:,4),cood(:,1):cood(:,1)+cood(:,3));
        %imag_face_s=225-imfill(225-imag_face_s,'holes');
        facebin=im2bw(imag_face_s,facethres);
    end
    figure(2);
    subplot(2,1,1);
    imshow(imag_face_s);
    subplot(2,1,2);
    imshow(facebin);
    [~,~,button]=ginput(1);
end

facebin=im2bw(imag_face_s,facethres);
figure(2);
imshow(facebin);

%select eye ROI manually
figure(1);
h2= impoly(gca);
wait(h2);
roi_eye=h2.createMask;
cood=round(h2.getPosition);
eyerect=[min(cood(:,2))-20,max(cood(:,2))+20,min(cood(:,1))-20,max(cood(:,1))+20];

% set binary threshold
% delete light spot
imag_eye_s=imag(eyerect(1):eyerect(2),eyerect(3):eyerect(4));
roi_eye_s=roi_eye(eyerect(1):eyerect(2),eyerect(3):eyerect(4));
inroi=imag_eye_s(roi_eye_s);
outroi=imag_eye_s(~roi_eye_s);
eyethres=mean(inroi)./mean(outroi);
eyethres=min(max(eyethres,0.01),0.99);
imagbin=im2bw(imag_eye_s,eyethres);
figure(3);
imshow(imagbin);
[~,~,button]=ginput(1);

while button~=49
    [~,~,button]=ginput(1);
    if button==43
        eyethres=min(eyethres+0.025,1);
        imagbin=im2bw(imag_eye_s,eyethres);
    elseif button==45
        eyethres=eyethres-0.025;
        imagbin=im2bw(imag_eye_s,eyethres);
    end
    figure(3);
    imshow(imagbin);
    [~,~,button]=ginput(1);
end

imagbin=im2bw(imag_eye_s,eyethres);
figure(3);
imshow(imagbin);


%now adjust the threshold for pupil detection
figure(4);
pupilthres=eyethres/2;
imag = rgb2gray(read(mov,2000));
imag_eye_s=imag(min(cood(:,2))-20:max(cood(:,2))+20,min(cood(:,1))-20:max(cood(:,1))+20);
imag_eye_s=255-imfill(255-imag_eye_s);
imag_eye_s = imgaussfilt(imag_eye_s,4);

pupilbin=im2bw(imag_eye_s,pupilthres);
imshow(pupilbin);
[~,~,button]=ginput(1);

while button~=49
    if button==43
        pupilthres=min(pupilthres+0.005,1);
        pupilbin=im2bw(imag_eye_s,pupilthres);
    elseif button==45
        pupilthres=pupilthres-0.005;
        pupilbin=im2bw(imag_eye_s,pupilthres);
    elseif button==50
        imag=rgb2gray(read(mov,round(2000+1000*rand())));
        imag_eye_s=imag(min(cood(:,2))-20:max(cood(:,2))+20,min(cood(:,1))-20:max(cood(:,1))+20);
        imag_eye_s=255-imfill(255-imag_eye_s);
        imag_eye_s= imgaussfilt(imag_eye_s,4);
        pupilbin=im2bw(imag_eye_s,pupilthres);
    end
    figure(4);
    imshow(pupilbin);
    [~,~,button]=ginput(1);
end

pupilbin=im2bw(imag_eye_s,pupilthres);
figure(4);
imshow(pupilbin);

close all;