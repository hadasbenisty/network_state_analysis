function newP = diffmap2clusters(matr)
if size(matr,3)>1
    ispos = squeeze(matr(77,182,:)>0);
    for ai = 1:size(matr,3)
        if ispos(ai)
            P=double(matr(:,:,ai));
        else
            P=double(-matr(:,:,ai));
        end
        P(isnan(matr(:,:,ai)))=nan;
        newP(:,:,ai) = P;
    end
else
    
    ispos = matr(1,:)>0;
    for ai = 1:size(matr,2)
        if ispos(ai)
            newP(:,ai)=double(matr(:,ai));
        else
            newP(:,ai)=double(-matr(:,ai));
        end
        
    end
end