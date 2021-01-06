function statenum = statename2num( statename)
switch statename
    case 'low_pup_q'
        statenum=1;
    case 'high_pup_q'
        statenum=2;
    case 'high_pup_l'
        statenum=3;
end