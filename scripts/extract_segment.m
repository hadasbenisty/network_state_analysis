function Y = extract_segment(t, X, segments)
Y=[];
L=min(length(t), size(X,2));
X=X(:,1:L);
t=t(1:L);
toremove = sum(segments>t(end) | segments<t(1),2);
segments=segments(toremove==0,:);
for i = 1:size(segments,1)
    if t(1)>segments(i, 1)
        continue;
    end
    if t(end) < segments(i, 2)
        continue;
    end
   ind1 = findClosestDouble(t, segments(i, 1));
   ind2 = findClosestDouble(t, segments(i, 2));
    Y = cat(2, Y, X(:, ind1:ind2));
    
end