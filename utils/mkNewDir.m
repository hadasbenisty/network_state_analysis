function mkNewDir(newdir)

if ~isempty(newdir) && ~exist(newdir, 'dir')
    mkdir(newdir);
end