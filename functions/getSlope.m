function slopeData = getSlope(X)

for pari = 1:size(X,1)
    p = polyfit(1:size(X,2),X(pari,:),1);
    slopeData(pari) = p(1); %#ok<AGROW>
end
end
