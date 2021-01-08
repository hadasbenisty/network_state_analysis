function animals_db = get_animals_meta_data_by_csv(csvfile)
if exist('../meta_data_processing/antara_exp/animals_db.mat','file')
    load('../meta_data_processing/antara_exp/animals_db.mat', 'animals_db');
else
    if nargin == 0
        csvfile = '../meta_data_processing/antara_exp/Processing_Pipeline_Full.csv';
    end
    T = readtable(csvfile);
    
    animal_lut = unique(T.animal);
    sessionsids_lut = unique(T.sessid);
    cohort_lut = unique(T.Cohort);
    type_lut = unique(T.type);
    N = length(T.animal);
    animal_list = zeros(N,1);
    sessionsid_list = zeros(N,1);
    cohort_list = zeros(N,1);
    folder_list = cell(N,1);
    type_list = zeros(N,1);
    
    isgoodpupil_list = ones(size(T.ispupilgood));
    isgoodpupil_list(T.ispupilgood>1) = 2;
    isgoodpupil_list(T.ispupilgood<=1) = 1;
    isgoodpupil_lut = {'GOOD','BAD'};
    
    isimagingood_list = ones(size(T.isimagingood));
    isimagingood_list(T.isimagingood>1) = 2;
    isimagingood_list(T.isimagingood<=1) = 1;
    isimagingood_lut = {'GOOD','BAD'};
    for k = 1:N
        animal_list(k) = find(strcmp(animal_lut, T.animal{k}));
        sessionsid_list(k) = find(strcmp(sessionsids_lut, T.sessid{k}));
        cohort_list(k) = find(strcmp(cohort_lut, T.Cohort{k}));
        type_list(k) = find(strcmp(type_lut, T.type{k}));
        folder_list{k} = T.directory{k};
    end
    
    animals_db.animal_list=animal_list;
    animals_db.sessionsid_list=sessionsid_list;
    animals_db.cohort_list=cohort_list;
    animals_db.folder_list=folder_list;
    animals_db.type_list=type_list;
    animals_db.animal_lut=animal_lut;
    animals_db.sessionsids_lut=sessionsids_lut;
    animals_db.cohort_lut=cohort_lut;
    animals_db.type_lut=type_lut;

animals_db.isgoodpupil_lut=isgoodpupil_lut;
animals_db.isgoodpupil_list=isgoodpupil_list;
    save('../meta_data_processing/antara_exp/animals_db.mat', 'animals_db');
end