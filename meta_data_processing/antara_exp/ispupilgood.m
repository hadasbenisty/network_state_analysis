function ispup_good =ispupilgood(animals_db, k)


ispup_good = animals_db.isgoodpupil_list(k)==find(strcmp(animals_db.isgoodpupil_lut, 'GOOD'));
   