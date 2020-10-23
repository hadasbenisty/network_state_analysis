addpath(genpath('../../../utils/'));
cd('X:\Lav\ProcessingDirectory')

addpath(genpath('X:\Lav\network_state_analysis\scripts\parcor_undirected'))
animals={'xt','xs','xu','xx','xz','xw'};

for ir=1:length(animals)
    demo_all_highfacelowface(char(animals(ir)), 'X:\Lav\ProcessingDirectory\', false);
end