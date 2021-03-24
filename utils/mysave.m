function mysave(f,filename, est)
if nargin == 2
    est = 'all';
end
[~,~,EXT] = fileparts(filename);
if ~isempty(EXT)
    est = EXT(2:end);
end
if isempty(f)
    f=gcf;
end
switch est
    case 'all'
        print( f, '-depsc', [filename '.eps']);
        saveas(f,[filename '.fig'],'fig');
        saveas(f,[filename '.tif'],'tif');
saveas(f,[filename '.jpg'],'jpg');
        % s=getframe(f);imwrite(s.cdata, filename,'Compression','none', 'Resolution',300 );
%         print( f, '-depsc', [filename '.eps']);
        % print(f,[filename '.pdf'],'-dpdf');
    case 'fig'
        saveas(f,[filename '.fig'],'fig');
    case 'tif'
        saveas(f,[filename '.tif'],'tif');
    case 'eps'
        print( f, '-depsc', [filename '.eps']);
end