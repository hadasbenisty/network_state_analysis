function plot_modeling_final
mkNewDir('x:\Hadas\Meso-imaging\GRABS_results\Figures\arousal_modeling');
 %plot_gfp;
% plot_2p_ca
plot_1p_ca;

end
function plot_gfp
figsfolder = 'x:\Hadas\Meso-imaging\GRABS_results\Figures\arousal_modeling';
proc_output_folder = 'X:\Hadas\Meso-imaging\GRABS_results\processing_dirgfp\';
inputFolder='x:\GRABS_Data\Analyzed_SVDMethodPatch14\GFP\';
addpath('../load_network_analysis_results/');

folders = dir(inputFolder);
folders=folders(3:end);
folders=folders([folders.isdir]);
folders=folders([folders.isdir]);
folders=folders(2:end);
winsizeSec = [0.5 1 3 5 10 15];
    [R2_phi,C_phi, R2_x,C_x,ttls]=loadRes_gfp_activity2arousal(proc_output_folder, folders, winsizeSec);

N=6;
J=20;
%% R
R2_phi=nanmean(R2_phi,5);
R2_x=nanmean(R2_x,4);
Nphi = sum(~isnan(squeeze(R2_phi(1,80,:,:))),2);
M_activity_phi = squeeze(nanmean(R2_phi(:,J,:,:),4));
S_activity_phi = squeeze(nanstd(R2_phi(:,J,:,:),[],4))/sqrt(N);
M_activity_x  =nanmean(R2_x,3);
S_activity_x  =nanstd(R2_x,[],3);

% R2_phi_reg = loadRes_gfp_network_reg(proc_output_folder, folders, winsizeSec, -1);
% Nphi2 = sum(~isnan(squeeze(R2_phi_reg(1,80,:,:))),2);
% M_net_phi = squeeze(nanmean(R2_phi_reg(:,J,:,:),4));
% S_net_phi = squeeze(nanstd(R2_phi_reg(:,J,:,:),[],4))/sqrt(N);


figure;
for a=1:3%length(ttls)
subplot(3,1,a);
M(:,1) = M_activity_x(a,:);
M(:,2) = M_activity_phi(a,:);
%M(:,3) = M_net_phi(a,:);
S(:,1) = S_activity_x(a,:);
S(:,2) = S_activity_phi(a,:);
%S(:,3) = S_net_phi(a,:);

barwitherr(S([3 5],:)',M([3 5],:)');
set(gca,'XTickLabel',{ 'x(t)','\phi_x','\phi_c'});
ylabel(ttls{a});
end
legend('3sec','10sec');

%saveas(gcf,fullfile(figsfolder,'arousal_by2p_bars'))

figure;plot(nanmean(R2_phi(:,:,end,:),4)')
legend(ttls);
xlabel('Embedded Components');ylabel('R^2');

end

function plot_2p_ca
figsfolder = 'x:\Hadas\Meso-imaging\GRABS_results\Figures\arousal_modeling';
proc_output_folder = 'X:\Hadas\Meso-imaging\GRABS_results\processing_dir2p\';
inputFolder='X:\Andrew\2Ponly\good spont sessions';
addpath('../load_network_analysis_results/');
overall_res_file = fullfile(proc_output_folder, 'activity', 'activity2arousal_2p.mat');
folders = dir(inputFolder);
folders=folders(3:end);
folders=folders([folders.isdir]);
winsizeSec = [0.5 1 3 5 10 15];
if isfile(overall_res_file)&&0
    load(overall_res_file, 'R2_phi',   'R2_x', 'ttls');
    
else
    [R2_phi,C_phi, R2_x,C_x,ttls]=loadRes_2p_activity2arousal(proc_output_folder, folders, winsizeSec);
    save(overall_res_file, 'R2_phi', 'C_phi',  'R2_x', 'C_x' ,'ttls');
end
N=6;
J=80;
%% R
Nphi = sum(~isnan(squeeze(R2_phi(1,80,:,:))),2);
M_activity_phi = squeeze(nanmean(R2_phi(:,J,:,:),4));
S_activity_phi = squeeze(nanstd(R2_phi(:,J,:,:),[],4))/sqrt(N);
M_activity_x  =nanmean(R2_x,3);
S_activity_x  =nanstd(R2_x,[],3);

R2_phi_reg = loadRes_2p_network_reg(proc_output_folder, folders, winsizeSec, -1);
Nphi2 = sum(~isnan(squeeze(R2_phi_reg(1,80,:,:))),2);
M_net_phi = squeeze(nanmean(R2_phi_reg(:,J,:,:),4));
S_net_phi = squeeze(nanstd(R2_phi_reg(:,J,:,:),[],4))/sqrt(N);


figure;
for a=1:3%length(ttls)
subplot(3,1,a);
M(:,1) = M_activity_x(a,:);
M(:,2) = M_activity_phi(a,:);
M(:,3) = M_net_phi(a,:);
S(:,1) = S_activity_x(a,:);
S(:,2) = S_activity_phi(a,:);
S(:,3) = S_net_phi(a,:);

barwitherr(S([3 5],:)',M([3 5],:)');
set(gca,'XTickLabel',{ 'x(t)','\phi_x','\phi_c'});
ylabel(ttls{a});
end
legend('3sec','10sec');

%saveas(gcf,fullfile(figsfolder,'arousal_by2p_bars'))

figure;plot(nanmean(R2_phi_reg(:,:,end,:),4)')
legend(ttls);
xlabel('Embedded Components');ylabel('R^2');

end
function plot_1p_ca
parcellation_methods = {'dec_16' 'dec_64' 'dec_svd' 'Allen','LSSC'};

MiceAnalyze=5:10;
inputFolder='X:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice';

proc_output_folder='x:\Hadas\Meso-imaging\GRABS_results\processing_dir\';

figsfolder = 'x:\Hadas\Meso-imaging\GRABS_results\Figures\arousal_modeling';
addpath(genpath('../functions/'));
addpath(genpath('../visualization/'))

addpath(genpath('../utils'))

winsizeSec = [0.5 1 3 5 10 15];

tauv =  -1 ;
Np = length(parcellation_methods);
for np=1:Np
    R2joint_phi{np} = nan(4, length(winsizeSec), length(MiceAnalyze),10);
    R2activity_phi{np} = nan(4, 99,length(winsizeSec), length(MiceAnalyze),10);
    R2_activity_x{np} = nan(4, length(winsizeSec), length(MiceAnalyze),10);
    
    R2net_phi_reg{np} = nan(4, 99,length(winsizeSec), length(MiceAnalyze),10);
    R2net_phi_noreg{np} = nan(4, 99,length(winsizeSec), length(MiceAnalyze),10);
    R2net_c{np} = nan(4, length(winsizeSec), length(MiceAnalyze),10);
    R2net_x{np} = nan(4, length(winsizeSec), length(MiceAnalyze),10);
end
for i=1:length(MiceAnalyze)
    animal = MiceAnalyze(i);
    mainDir1 =fullfile(inputFolder, sprintf("grabAM%02d",(animal)), 'imaging with 575 excitation');
    folders = dir(mainDir1);
    folders=folders(3:end);
    for folder =1:length(folders)
        
        currfolder = fullfile(mainDir1, folders(folder).name);
        if  contains(currfolder, 'vis') ||contains(currfolder, 'PostDrug')%
            continue;
        end
         %% joint 
        for np=1:Np
            parcellation = parcellation_methods{np};
            for wi = 1:length(winsizeSec)
                switch parcellation
                    case 'Allen'
                        resfile=fullfile(proc_output_folder, 'joint', [ 'reg_-1_' parcellation '_green_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                    otherwise
                        resfile=fullfile(proc_output_folder, 'joint',  ['reg_-1_' parcellation '_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                end
                if ~isfile(resfile)
                    
                    continue;
                end
                try
                    res =load(resfile);
                    
                catch
                    disp(['problem loading ' resfile])
                    continue;
                end
             
                R2joint_phi{np}(:,  wi,i,folder) = max(0,res.R_phi);                
            end
            
           
        end
        %% activity no reg R2activity_phi
        for np=1:Np
            parcellation = parcellation_methods{np};
            for wi = 1:length(winsizeSec)
                switch parcellation
                    case 'Allen'
                        resfile=fullfile(proc_output_folder, 'activity', [  parcellation '_green_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                    otherwise
                        resfile=fullfile(proc_output_folder,  'activity', [parcellation '_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                end
                if ~isfile(resfile)
                    disp(resfile);
                    continue;
                end
                try
                    res =load(resfile);
                    if size(res.R_phi,1)==3
                        R2activity_phi{np}(1:3, :, wi,i,folder) = max(0,res.R_phi);
                        R2_activity_x{np}(1:3, wi,i,folder) = max(0,res.R_x);
                        continue;
                    end
                catch
                    disp(['problem loading ' resfile])
                    continue;
                end
                
                R2activity_phi{np}(:, :, wi,i,folder) = max(0,res.R_phi);
                R2_activity_x{np}(:, wi,i,folder) = max(0,res.R_x);
            end
            
            
        end
        
        %% net no reg R2net_phi_noreg
        for np=1:Np
            parcellation = parcellation_methods{np};
            for wi = 1:length(winsizeSec)
                switch parcellation
                    case 'Allen'
                        resfile=fullfile(proc_output_folder, 'network_analysis', [ 'reg_0_' parcellation '_green_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                    otherwise
                        resfile=fullfile(proc_output_folder, 'network_analysis', ['reg_0_' parcellation '_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                end
                if ~isfile(resfile)
                    
                    continue;
                end
                try
                    res =load(resfile);
                    if size(res.R_phi,1)==3
                        R2net_phi_noreg{np}(1:3, :, wi,i,folder) = max(0,res.R_phi);
                        R2net_c{np}(1:3, wi,i,folder) = max(0,res.R_c);
                        R2net_x{np}(1:3,  wi,i,folder) = max(0,res.R_mx);
                        continue;
                    end
                catch
                    disp(['problem loading ' resfile])
                    continue;
                end
                if ~isempty(res.R_c)
                    R2net_c{np}(:, wi,i,folder) = max(0,res.R_c);
                    R2net_x{np}(:,  wi,i,folder) = max(0,res.R_mx);
                end
                R2net_phi_noreg{np}(:, :, wi,i,folder) = max(0,res.R_phi);
            end
            
            taui=1;
            for wi = 1:length(winsizeSec)
                switch parcellation
                    case 'Allen'
                        resfile=fullfile(proc_output_folder, 'network_analysis', ['reg_', num2str(tauv(taui)) '_' parcellation '_green_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                    otherwise
                        resfile=fullfile(proc_output_folder,  'network_analysis', ['reg_', num2str(tauv(taui)) '_' parcellation '_'  folders(folder).name '_sims' num2str(winsizeSec(wi))  '.mat']);
                end
                if ~isfile(resfile)
                    %   disp(resfile)
                    continue;
                end
                
                try
                    res =load(resfile,'R_phi','C_phi');
                    if size(res.R_phi,1)==3
                        R2net_phi_reg{np}(1:3, :, wi,i,folder) = max(0,res.R_phi);
                        continue;
                    end
                catch
                    disp(['problem loading ' resfile])
                    continue;
                end
                
                R2net_phi_reg{np}(:, :, wi,i,folder) = max(0,res.R_phi);
            end
        end
    end
    
end
ttls = {'Pupil', 'Face Map (1)', 'Wheel','pca'};
Jc=20;
for np=1:Np
    contribution_correlation_c{np} = squeeze(R2net_phi_reg{np}(:,Jc,:,:,:))./R2joint_phi{np};
    
    contribution_correlation_c{np}=nanmean(contribution_correlation_c{np},4);
    contribution_correlation_a{np} = squeeze(R2activity_phi{np}(:,Jc,:,:,:))./R2joint_phi{np};
   
     
    contribution_correlation_a{np}=nanmean(contribution_correlation_a{np},4);
    
    R2joint_phi{np} = nanmean(R2joint_phi{np},4);
    R2net_phi_reg{np} = nanmean(R2net_phi_reg{np},5);
    R2net_phi_noreg{np} = nanmean(R2net_phi_noreg{np},5);
    R2net_c{np} = nanmean(R2net_c{np},4);
    R2net_x{np} = nanmean(R2net_x{np},4);
    R2activity_phi{np} = nanmean(R2activity_phi{np},5);
    R2_activity_x{np} = nanmean(R2_activity_x{np},4);
    N_phi(:,np) = sum(~isnan(squeeze(R2net_phi_noreg{np}(1,1,:,:))),2); %#ok<NASGU>
    N_reg(:,np) = sum(~isnan(squeeze(R2net_phi_reg{np}(1,1,:,:))),2);%#ok<NASGU>
    N_j(:,np) = sum(~isnan(squeeze(R2joint_phi{np}(1,:,:))),2);%#ok<NASGU>

    N_netphi(:,np) = sum(~isnan(squeeze(R2net_phi_noreg{np}(1,1,:,:))),2); %#ok<NASGU>
    N_netreg(:,np) = sum(~isnan(squeeze(R2net_phi_reg{np}(1,1,:,:))),2);%#ok<NASGU>
    N_net_c(:,np) = sum(~isnan(squeeze(R2net_c{np}(1,:,:))),2); %#ok<NASGU>
    N_activity_phi(:,np) = sum(~isnan(squeeze(R2activity_phi{np}(1,1,:,:))),2);%#ok<NASGU>
    
end


N=length(MiceAnalyze);
J=80;
l=1;ttl=[];Mcontribution=[];ttlc=[];lc=1;
for np=1:Np
    p=parcellation_methods{np};
    p(p=='_') = ' ';
    if contains(p, 'dec')
        p=p(5:end);
    end
    M(:,:,l) = nanmean(R2_activity_x{np},3);
    S(:,:,l) = nanstd(R2_activity_x{np},[],3)/sqrt(N);
    ttl{l} = ['activity ' p];
    l=l+1;
    
    
    M(:,:,l) = squeeze(nanmean(R2activity_phi{np}(:,J,:,:),4));
    S(:,:,l) = squeeze(nanstd(R2activity_phi{np}(:,J,:,:),[],4)/sqrt(N));
    ttl{l} = ['\phi_a ' p];
    l=l+1;
end
J=20;
for np=1:Np
    p=parcellation_methods{np};
    p(p=='_') = ' ';
    M(:,:,l) = squeeze(nanmean(R2net_c{np},3));
    S(:,:,l) = squeeze(nanstd(R2net_c{np},[],3)/sqrt(N));
    ttl{l} = ['corr ' p];
    l=l+1;
    
    M(:,:,l) = squeeze(nanmean(R2net_x{np},3));
    S(:,:,l) = squeeze(nanstd(R2net_x{np},[],3)/sqrt(N));
    ttl{l} = ['proj ' p];
    l=l+1;
    
    M(:,:,l) = squeeze(nanmean(R2net_phi_noreg{np}(:,J,:,:),4));
    S(:,:,l) = squeeze(nanstd(R2net_phi_noreg{np}(:,J,:,:),[],4)/sqrt(N));
    ttl{l} = ['\phi_c ' p];
    l=l+1;
    
    M(:,:,l) = squeeze(nanmean(R2net_phi_reg{np}(:,J,:,:),4));
    S(:,:,l) = squeeze(nanstd(R2net_phi_reg{np}(:,J,:,:),[],4)/sqrt(N));
    ttl{l} = ['\phi_c reg ' p];
    l=l+1;
     M(:,:,l) = squeeze(nanmean(R2joint_phi{np}(:,:,:),3));
    S(:,:,l) = squeeze(nanstd(R2joint_phi{np}(:,:,:),[],3)/sqrt(N));
    ttl{l} = ['\phi_j reg ' p];
    l=l+1;
    
    Mcontribution(:,:,lc) = nanmean(contribution_correlation_a{np},3);
    Scontribution(:,:,lc) = nanstd(contribution_correlation_a{np},[],3);
    ttlcont{lc} = ['\phi_a ' p];
    lc=lc+1;
    
    Mcontribution(:,:,lc) = nanmean(contribution_correlation_c{np},3);
    Scontribution(:,:,lc) = nanstd(contribution_correlation_c{np},[],3);
    ttlcont{lc} = ['\phi_c ' p];
    lc=lc+1;
end
%% contribution
figure;
for a=1:3
    subplot(3,1,a);
    inds=any(~isnan(squeeze(Mcontribution(a,:,:))));
    barwitherr(squeeze(Scontribution(a,:,inds)), squeeze(Mcontribution(a,:,inds)));
    %
    ylim([0 1]);set(gca,'XTickLabel', []);xlim([2.5 6.5])
    ylabel(ttls{a});xlabel('Window Size [sec]');
end
set(gca,'XTickLabel', winsizeSec);legend(ttlcont(inds),'Location','Best');
for np=1:np
alpha_v{np}=contribution_correlation_c{np}./(contribution_correlation_c{np}+contribution_correlation_a{np});
Ma(:,:,np) = nanmean(alpha_v{np},3);
Sa(:,:,np) = nanstd(alpha_v{np},[],3);
end

figure;
for a=1:4
    subplot(4,1,a);
    
    inds=any(~isnan(squeeze(Ma(a,:,:))));
    barwitherr(squeeze(Sa(a,:,inds)), squeeze(Ma(a,:,inds)));
    %
    ylim([0 1]);set(gca,'XTickLabel', [])
    ylabel(ttls{a});
end
set(gca,'XTickLabel',winsizeSec);xlabel('Window Size [sec]');
suptitle('Relative contribution of correlation');

figure;
clr = 'kgbr';
for a=1:4  
    shadedErrorBar(winsizeSec, squeeze(nanmean(R2net_phi_noreg{4}(a,20,:,:),4)),...
    squeeze(nanstd(R2net_phi_noreg{4}(a,20,:,:),[],4))/sqrt(N),'lineprops',clr(a));hold all
end
set(gca,'XTick',winsizeSec);
xlabel('window size [sec]');c=get(gca,'Children');
legend(c(end:-4:1),ttls);axis tight;ylabel('R^2');
figure;plot(nanmean(R2net_phi_noreg{4}(:,:,end,:),4)');
legend(ttls);
xlabel('Embedded Components');ylabel('R^2');
figure;
for a=1:4
    subplot(4,1,a);
    inds=any(~isnan(squeeze(M(a,:,:))));
    barwitherr(squeeze(S(a,:,inds)), squeeze(M(a,:,inds)));
    %
    ylim([0 1]);set(gca,'XTickLabel', []);
    ylabel(ttls{a});xlabel('Window Size [sec]');
end
set(gca,'XTickLabel', winsizeSec);legend(ttl(inds),'Location','Best');
figure;
for a=1:4
    subplot(4,1,a);
    inds=any(~isnan(squeeze(M(a,:,:))));
    b=barwitherr(squeeze(S(a,winsizeSec==15,inds)), squeeze(M(a,winsizeSec==15,inds)));
    set(gca,'XTick',1:sum(inds));
    set(gca,'XTickLabel', ttl(inds));ylim([0 .9]);
    title(ttls{a});ylabel('R^2');
end
inds=find(inds);
inds=inds(4:end);
figure;
b=barwitherr(squeeze(S(:,winsizeSec==10,inds))', squeeze(M(:,winsizeSec==10,inds))');
set(gca,'XTick',1:sum(inds));
set(gca,'XTickLabel', ttl(inds));ylim([0 .9]);
legend(ttls);ylabel('R^2');
title('Window Size 10secs');
return;
figure;
for a=1:4
    M(:,1) = squeeze(R2green_phiM(a,20,:));
    M(:,2) = squeeze(R2green_phi_regM(a,20,:));
    
    S(:,1) = squeeze(R2green_phiS(a,20,:));
    S(:,2) = squeeze(R2green_phi_regS(a,20,:));
    %figure;plot(winsizeSec, M);legend({'0','0.01','0.1','1','2','4'});
    
    subplot(4,1,a);
    shadedErrorBar(winsizeSec, M(:,1),S(:,1));hold all;
    %shadedErrorBar(winsizeSec, M(:,2),S(:,2),'lineProps','b');
    %shadedErrorBar(winsizeSec, M(:,3),S(:,3),'lineProps','c');
    shadedErrorBar(winsizeSec, M(:,2),S(:,2),'lineProps','g');
    %shadedErrorBar(winsizeSec, M(:,5),S(:,5),'lineProps','y');
    %shadedErrorBar(winsizeSec, M(:,6),S(:,6),'lineProps','m');
    %shadedErrorBar(winsizeSec, M(:,7),S(:,7),'lineProps','r');
    %shadedErrorBar(winsizeSec, M(:,8),S(:,8),'lineProps','c');
    c=get(gca,'Children');
    title(ttls{a});ylabel('R^2');xlabel('Window Size [sec]');
end
legend(c(end:-4:1),{'no reg','\tau-svd'})%'interp',

%legend(c(end:-4:1),{'no reg','pca','svd','\tau=1e-4','\tau=0.01','\tau=0.1','\tau=1'})%'interp',
figure;
for a=1:4
    M(:,1) = squeeze(R2green_phiM(a,20,:));
    M(:,2) = squeeze(R2green_phi_regM(a,20,:));
    
    S(:,1) = squeeze(R2green_phiS(a,20,:));
    S(:,2) = squeeze(R2green_phi_regS(a,20,:));
    %figure;plot(winsizeSec, M);legend({'0','0.01','0.1','1','2','4'});
    
    subplot(4,1,a);
    barwitherr(S, M);hold all;
    legend({'no reg','\tau-svd'},'Location','Best');
    set(gca,'XTickLabel', winsizeSec);ylim([0 1]);
    title(ttls{a});ylabel('R^2');xlabel('Window Size [sec]');
end
if strcmp(parcellation, 'Allen')
    uiopen('/gpfs/gibbs/pi/cardin-higley/Hadas/results/figs/network/R2comparison.fig',1);
    c=get(gcf,'Children');
    axes(c(3));
    hold all;
    shadedErrorBar(winsizeSec, squeeze(R2green_phi_regM(1,20,:)),squeeze(R2green_phi_regS(1,20,:)),'lineProps','g');
    axes(c(2));
    hold all;
    shadedErrorBar(winsizeSec, squeeze(R2green_phi_regM(2,20,:)),squeeze(R2green_phi_regS(2,20,:)),'lineProps','g');
    axes(c(1));
    hold all;
    shadedErrorBar(winsizeSec, squeeze(R2green_phi_regM(3,20,:)),squeeze(R2green_phi_regS(3,20,:)),'lineProps','g');
    set(gca,'XTickLabel', winsizeSec);ylim([0 1]);
    c=get(gca,'Children');
    legend(c(end:-4:1),{'activity','correlation'},'Location','Best');
    ylabel('R^2');xlabel('Window Size [sec]');
    saveas(gcf,fullfile(figsfolder, [parcellation 'R2comparison1.fig']));
end
end


