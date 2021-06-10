function pupilframe = extract_pupil_from_data(pupilframe)



% decide high/low polarity
if nanmean(pupilframe)>1.5 % low polarity
    pupilframe=-(pupilframe-2.9);
end
[pupon,pupoff]=squaredetect(pupilframe,0.05*max(pupilframe)); %max value=2.9
pupon=pupon(1:length(pupoff));
if pupon(100)-pupon(99)>5000/30*1.7
    maxi=max(pupilframe); %step amp of each frame
    for i=1:length(pupon)-1
        if round((pupoff(i+1)+pupoff(i))/2)<pupon(end)
            pupilframe(round((pupon(i+1)+pupon(i))/2):round((pupoff(i+1)+pupoff(i))/2))=maxi;
        end
    end
end
pupilframe=(pupilframe-nanmin(pupilframe)>0.5);