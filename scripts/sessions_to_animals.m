function [animal_values,animaltype_list,animal_idlist] = sessions_to_animals(session)

if ndims(session)==2&&((size(session,1)==1&&size(session,2)==85)||(size(session,2)==1&&size(session,1)==85))
    animals_db = get_animals_meta_data_by_csv;
    %%sessions to animals
    counter=1; %counter for index to fill growing dataset
    animaltype_list=[]; %vector of animal type labels ordered 
    animal_values=[]; %vector of animal values as per the function input ordered
    animal_idlist=[]; %vector of animal IDS ordered
    %animal_list=unique(animals_db.animal_list(animals_db.toinclude_list==3)); %allunique animals
 %   animal_list=[1 2 3 4 5 6 10 11 12 13 14 19 21 22 24 25 27 28 33 35 36 37];
    animal_list=[1 2 4 6 10 11 12 13 14 19 21 22 24 25 27 28 33 35 36 37];

    for j=1:length(animal_list) %iterate through all unique animals
        curr_animalid=animal_list(j);  %unique animal ID for each iteration
        if curr_animalid==8||curr_animalid==9||curr_animalid==23 %these animals only have 1 state
            continue;
        end
        allsessions=session(find(animals_db.animal_list==curr_animalid)); %find sessions per animal ID
        if isempty(allsessions)||all(isnan(allsessions(:))) %if no valid sessions, skip animal ID
            %continue;
            allsessions=NaN;
        end
        animal_values(counter)=nanmean(allsessions); %average across sessions per animal
        animaltype=animals_db.type_list(animals_db.animal_list==curr_animalid); %get animal type for animalID
        animaltype_list(counter)=animaltype(1);
        animal_idlist(counter)=curr_animalid; %store animal ID in same order as type and values
        counter=counter+1;
        clearvars allsessions
    end
else
   % disp('not right length')
end

%%matrix form
if ndims(session)==3
    disp('larger than vector')
%     animals_db = get_animals_meta_data_by_csv;
%     %%sessions to animals
%     counter=1; %counter for index to fill growing dataset
%     animaltype_list=[]; %vector of animal type labels ordered 
%     animal_values=[]; %vector of animal values as per the function input ordered
%     animal_idlist=[]; %vector of animal IDS ordered
%     animal_list=unique(animals_db.animal_list); %allunique animals
%     for j=1:length(unique(animals_db.animal_list)) %iterate through all unique animals
%         curr_animalid=animal_list(j);  %unique animal ID for each iteration
%         if curr_animalid==8||curr_animalid==9||curr_animalid==23 %these animals only have 1 state
%             continue;
%         end
%         allsessions=session(find(animals_db.animal_list==curr_animalid)); %find sessions per animal ID
%         if isempty(allsessions)||all(isnan(allsessions(:))) %if no valid sessions, skip animal ID
%             continue;
%         end
%         animal_values(counter)=nanmean(allsessions); %average across sessions per animal
%         animaltype=animals_db.type_list(animals_db.animal_list==curr_animalid); %get animal type for animalID
%         animaltype_list(counter)=animaltype(1);
%         animal_idlist(counter)=curr_animalid; %store animal ID in same order as type and values
%         counter=counter+1;
%         clearvars allsessions
%     end
else
   % disp('not right length')
end

end

