
function makeMovie(videoName,finvideo,climit,frameRate)
%set up the movie 
writerObj = VideoWriter(videoName); 
writerObj.FrameRate = 10;
open(writerObj); figure;
for i=1:size(finvideo,3)
imagesc(finvideo(:,:,i));caxis(climit); colormap('jet'); colorbar; 
pause(1/frameRate)
frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
writeVideo(writerObj, frame)
end 
close(writerObj);
end 


