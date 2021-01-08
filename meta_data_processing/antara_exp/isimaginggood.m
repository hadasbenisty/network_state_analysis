function isimaging_good =isimaginggood(animals_db, k)
isimaging_good = animals_db.isimagingood_list(k)==find(strcmp(animals_db.isimagingood_lut, 'GOOD'));
