function plot_cr_Res(phig,cg,mxg,phib,cb,mxb, N,winsizeSec,ttls, clm)



Mphig = nanmean(phig, 3);
SEMphig = nanstd(phig, [],3)./sqrt(N);

Mcg = nanmean(cg, 3);
SEMcg = nanstd(cg, [],3)./sqrt(N);

Mmxg = nanmean(mxg, 3);
SEMmxg = nanstd(mxg, [],3)./sqrt(N);

Mphib = nanmean(phib, 3);
SEMphib = nanstd(phib, [],3)./sqrt(N);

Mcb = nanmean(cb, 3);
SEMcb = nanstd(cb, [],3)./sqrt(N);

Mmxb = nanmean(mxb, 3);
SEMmxb = nanstd(mxb, [],3)./sqrt(N);


for i=1:3
    subplot(2,3,i);
    M = [ Mcg(i,:); Mmxg(i,:);Mphig(i,:) ];
    SEM = [ SEMcg(i,:); SEMmxg(i,:);SEMphig(i,:)];
barwitherr(SEM',M');
set(gca,'XTickLabel', winsizeSec);title(['Green ' ttls{i}]);
xlabel('Window [sec]');ylim(clm);

subplot(2,3,i+3);
    M = [ Mcb(i,:); Mmxb(i,:);Mphib(i,:) ];
    SEM = [ SEMcb(i,:); SEMmxb(i,:);SEMphib(i,:)];
barwitherr(SEM',M');
set(gca,'XTickLabel', winsizeSec);title(['Blue ' ttls{i}]);
xlabel('Window [sec]');ylim(clm);
end
legend({'C(t)','proj{C(t)}','\phi(t)'}, 'Location','Best');