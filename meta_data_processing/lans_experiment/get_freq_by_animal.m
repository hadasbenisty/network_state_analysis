function fsimaging = get_freq_by_animal(animalName)

switch animalName
        case {'xu','xv','xt','xs'}
            fsimaging=33;
        otherwise
            fsimaging=10;
end