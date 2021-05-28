function [h,p]=corr_perm_test(mat1, mat2, mat1sh, mat2sh, inds)


mat1M = nanmean(mat1(:,:,inds),3);
mat2M = nanmean(mat2(:,:,inds),3);
diffM = mat1M-mat2M; %take the difference of correlation matrix mean across animals
mat1M_sh = nanmean(mat1sh(:,:,:,inds),4);
mat2M_sh = nanmean(mat2sh(:,:,:,inds),4);
diffM_sh = mat1M_sh-mat2M_sh; %take the differe
pval = mean(abs(diffM_sh)>=abs(diffM),3);

pval=pval-diag(diag(pval));
pval(pval==0)=1/1000;
[h,~,~,p]=fdr_bh(pval,0.01,'dep'); %benjamini
% x = reshape(diffM_sh,23*23,[]);
% y = reshape(diffM,23*23,[]);
% pp=ttest(bsxfun(@minus,x,y)',0.05,'tail','left');
% figure;imagesc(reshape(pp,23,23))
% ttest(mean(abs(diffM_sh)>=abs(diffM),3);

