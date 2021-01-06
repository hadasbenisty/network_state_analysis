function [Allparcells, parcels, maskAllen, maskAllenRegion, ROI_list] = getParcellsByLansAllansAtlas
ROI_list=[];
parcellsMat = 'AllenParcels.mat';
parcellsMat='X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat';
if exist(parcellsMat, 'file')
%     load(parcellsMat,  'Allparcells', 'parcels');
    load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat', 'parcells_new');
    parcels=parcells_new;
else
    load('allregions.mat', 'allregions');
    [~,TXT] = xlsread('../parcellation/AllenParcellationLan/allanParcellationTiffs/subregion_list.csv');
    N_parcells = length(TXT)-1;
    Allparcells = zeros(size(allregions));
    indicators = zeros([size(allregions) N_parcells]);
    names = cell(N_parcells, 1);
    description= cell(N_parcells, 1);
    files = dir('../AllenParcellationLan/allanParcellationTiffs/*.tif');
    for n = 1:length(files)
        a = imread(fullfile(files(n).folder, files(n).name));
        a = imresize(a, size(allregions));
        ind = find(strcmpi(TXT(:,1), files(n).name(1:end-4)));
        indicators(:,:,ind-1) = double(rgb2gray(a)>0);
        names{ind-1} = TXT{ind,1};
        description{ind-1} = TXT{ind,2};
        Allparcells = Allparcells + double(rgb2gray(a)>0)*n;
    end
    parcels.names = names;
    parcels.description = description;
    parcels.indicators = indicators;
    save(parcellsMat,  'Allparcells', 'parcels');
end
Allparcells = parcels.CombinedParcells;

if nargout==0
    figure;imagesc(Allparcells);
    
    placeRoiLabels(parcels);
end
parcels.regionNum = zeros(length(parcels.names),1);
parcels.regionNum(1:16) = 1; % visual
parcels.regionNum(17:22) = 5;% retrosplenial
parcels.regionNum(27:28) = 3; % temporal
parcels.regionNum(29:30) = 4;% auditory
parcels.regionNum(31:32) = 2; % parietal
parcels.regionNum(33:48) = 6; % somatosensory
parcels.regionNum(49:52) = 7; % somatomotor
parcels.regionNum(53:56) = 0; % rest



for k=1:length(parcels.regionNum)
    switch parcels.regionNum(k)
         case 0
            parcels.regionStr{k} = 'rest';
        case 1
            parcels.regionStr{k} = 'Visual';
        case 2
            parcels.regionStr{k} = 'Par';
        case 3
            parcels.regionStr{k} = 'Temp';
        case 4
            parcels.regionStr{k} = 'Aud';
        case 5
            parcels.regionStr{k} = 'RS';
        case 6
            parcels.regionStr{k} = 'S-S';
        case 7
            parcels.regionStr{k} = 'Motor';
    end
end
[~,maskAllen] = getAllenClusteringLabels([], [], parcels);


maskAllenRegion=zeros(256);
for pi=1:56
    maskAllenRegion(maskAllen==pi) = parcels.regionNum(pi);
end
if nargin == 5
nROI=56;
ss=size(maskAllenRegion);
ROI_list=[];
for i=1:nROI
    tempROI.pixel_list=find(parcels.indicators(:,:,i)==1);
    ROImap = parcels.indicators(:,:,i);
    tempROI.boundary_list=find(bwperim(ROImap,4));
    
    
    
    [I1,J1]=ind2sub(ss,tempROI.pixel_list);
    tempROI.centerPos=[mean(I1),mean(J1)];
    
    
    
    tempROI.seedPos=tempROI.centerPos;
    
    
    
    tempROI.name=parcels.names{i};
    
    
    ROImap_fill=imfill(ROImap,'holes');
    
    tempROI.fill_list=find(ROImap_fill);
    
    
    
    tempROI.fmean=0;
    tempROI.fmean_fill=0;
    ROI_list=[ROI_list,tempROI];

end
end


L=zeros(13+256,256,56);

    L(1:256,:,:) = parcels.indicators;
parcels.indicators=L(14:end,1:256,:);
Allparcells = zeros(size(size(parcels.indicators,1)));

for n = 1:size(parcels.indicators,3)   
    Allparcells = Allparcells + double(parcels.indicators(:,:,n))*n;
end
