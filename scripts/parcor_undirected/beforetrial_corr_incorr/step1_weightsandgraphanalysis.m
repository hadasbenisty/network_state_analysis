addpath(genpath('../../../utils/'));
cd('X:\Lav\ProcessingDirectory');

addpath(genpath('X:\Lav\network_state_analysis\scripts\parcor_undirected'))
animals={'xt','xs','xu','xx','xz','xw'};

for ir=1:length(animals)
    demo_all_beforetrial_corrincorr('X:\Lav\ProcessingDirectory_spet29\', char(animals(ir)), true);
    demo_all_beforetrial_corrincorr('X:\Lav\ProcessingDirectory\', char(animals(ir)), false);
end