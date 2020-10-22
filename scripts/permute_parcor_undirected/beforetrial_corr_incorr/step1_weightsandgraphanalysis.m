addpath(genpath('../../../utils/'));
cd('X:\Lav\ProcessingDirectory')

addpath(genpath('X:\Lav\network_state_analysis\scripts\permute_parcor_undirected'))
animals={'xt','xu','xx','xz','xw','xs'};
for ir=1:length(animals)
    perm_beforetrial_corincorr(char(animals(ir)))
    close all
end