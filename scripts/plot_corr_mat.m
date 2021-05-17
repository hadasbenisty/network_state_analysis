function plot_corr_mat(C,parcellnames,clm)
imagesc(tril(C),clm);
% axis equal
c=gca;
c.Box='off';
if isempty(parcellnames)
c.XAxis.Visible='off';
c.YAxis.Visible='off';
else
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);
end
ylim(c.XAxis.Limits);

