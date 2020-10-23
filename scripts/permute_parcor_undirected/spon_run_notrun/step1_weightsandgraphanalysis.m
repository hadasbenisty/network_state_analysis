addpath(genpath('../../../utils/'));
cd('X:\Lav\ProcessingDirectory')

addpath(genpath('X:\Lav\network_state_analysis\scripts\permute_parcor_undirected'))
animals={'xt','xs','xu','xx','xz','xw'};

for ir=1:length(animals)
    perm_runnotrun(char(animals(ir)));
end