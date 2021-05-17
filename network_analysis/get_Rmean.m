function [mRiemannianMean,tC2] = get_Rmean(weights)
Np = sqrt(size(weights,2));
tC2=zeros(Np, Np, size(weights,3));
for tt=1:size(weights,3)
    tC2(:,:,tt)=reshape(weights(1,:,tt),[Np Np])+eye(Np);
end



mRiemannianMean = RiemannianMean(tC2);
