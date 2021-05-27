function [animal_values,animal_list] = sessions_to_animals_3states(session_3states)
if size(session_3states,1)*size(session_3states,2)==225
    animals_db = get_animals_meta_data_by_csv;
    %%sessions to animals
    counter=1;
    animal_list=[];
    animal_values=[];
    allstates=[3,2,1];
    
    for state_i=1:3
        state=allstates(state_i);
        
        for j=1:length(unique(animals_db.animal_list))
            allsessions=session_3states(find(animals_db.animal_list==j));
            if isempty(allsessions)|all(isnan(allsessions(:)))
                continue;
            end
            animal_values(counter)=nanmean(allsessions);
            animaltype=animals_db.type_list(animals_db.animal_list==j);
            animal_list(counter)=animaltype(1);
            counter=counter+1;
            clearvars allsessions
        end
else
    disp('not right length')
end
end

