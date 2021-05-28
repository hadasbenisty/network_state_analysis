function [animal_values,animal_idlist] = sessions_to_animals_type(session,type)

if ndims(session)==2
    animals_db = get_animals_meta_data_by_csv;
    %%sessions to animals
    counter=1; %counter for index to fill growing dataset
    animal_values=[]; %vector of animal values as per the function input ordered
    animal_idlist=[]; %vector of animal IDS ordered
    animal_list=unique(animals_db.animal_list(animals_db.type_list==type&animals_db.toinclude_list==3)); %allunique animals
    
    new_folder_list=animals_db.folder_list(animals_db.type_list==type&animals_db.toinclude_list==3);
    new_animal_list=animals_db.animal_list(animals_db.type_list==type&animals_db.toinclude_list==3);
    for j=1:length(animal_list) %iterate through all unique animals
        curr_animalid=animal_list(j);  %unique animal ID for each iteration
        if curr_animalid==8||curr_animalid==9||curr_animalid==23 %these animals only have 1 state
            continue;
        end
        %need to order sessions by animals of this type that are 
        allsessions=session(find(new_animal_list==curr_animalid)); %find sessions per animal ID
        if isempty(allsessions)||all(isnan(allsessions(:))) %if no valid sessions, skip animal ID
            continue;
        end
        animal_values(counter)=nanmean(allsessions); %average across sessions per animal
        animal_idlist(counter)=curr_animalid; %store animal ID in same order as type and values
        counter=counter+1;
        clearvars allsessions
    end
else
   % disp('not right length')
end


end

