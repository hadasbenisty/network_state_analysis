function [n,dayss,lengthdays,midpoint] = animaltodays(animal)
if strcmp(animal,'xs')
    n=strcat(animal,num2str(32));
    dayss=[15 16 17 18 19 20 21 24 25 27 29 30 31 32];% 18
elseif strcmp(animal,'xv')
    n=strcat(animal,num2str(34));
    dayss=[23 24 25 26 27 28 29 31 32 33 34];
elseif strcmp(animal,'xu')
    n=strcat(animal,num2str(37));
    dayss=[20 22 24 26 27 29 30 32 33 34 35 36 37];% add 31
elseif strcmp(animal,'xx')
    n=strcat(animal,num2str(35));
    dayss=[12 14 15 16 17 18 19 20 21 22 23 25 26 27 29 30 33 34 35];
elseif strcmp(animal,'xz')
    n=strcat(animal,num2str(36));
    dayss=[18 19 20 21 22 23 24 25 26 27 28 29 30 33 34 35 36];% remove 19 25 27 
elseif strcmp(animal,'xt')
    n=strcat(animal,num2str(32));
    dayss=[15 17 18 19 21 22 23 25 26 27 28 29 30 32];% remove  25 27
elseif strcmp(animal,'xw')
    n=strcat(animal,num2str(36));
    dayss=[13 14 15 16 17 18 19 21 22 23 24 25 26   27 28 29 30 33 34 35 36];% remove 16 21 22 23 24  25 26 35 36 replace 13 with 12
end
lengthdays=length(dayss);
midpoint=round(lengthdays/2);
end

