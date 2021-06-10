
% say you have a blue movie and a uv movie (you don't really need the uv..)
[R, C, T] = size(blue_mov);
% transfrom
[tform, blue_mov_trans, uv_mov_trans] = parcellateData(blue_mov, uv_mov);
% get parcellation
[Allparcells, parcells] = getParcellsByLansAllansAtlas;
% reshape to pixels over time
blue_mov_trans_ii = reshape(blue_mov_trans, R*C , T);
uv_mov_trans_ii = reshape(uv_mov_trans, R*C , T);
blue_parcells = zeros(R*C, T);
% get parcelled data
for parcel_i = 1:length(parcells.names)
    K = find(parcells.indicators(:,:,parcel_i)==1);     
    blue_parcells(parcel_i, :) = nanmean(blue_mov_trans_ii(K,:),1);   
    uv_parcells(parcel_i, :) = nanmean(uv_mov_trans_ii(K,:),1);   
end

for parcel_i = 1:length(parcells.names)
    K = find(parcells.indicators(:,:,parcel_i)==1);     
    dF_parcells(parcel_i, :) = nanmean(X(K,:),1);   
end
plot(
t_imaging = 33 Hz 
figure;plot(dF_parcells(1,:));