function test4Bugs

condition1='run';condition2='notrun';
findIdentical(condition1, condition2)
condition1='pupilhigh';
condition2='pupillow';
findIdentical(condition1, condition2)

condition1='trial_notrun_correct';
condition2='trial_notrun_incorrect';
findIdentical(condition1, condition2)

condition1='trial_lowpup_correct';
condition2='trial_lowpup_incorrect';
findIdentical(condition1, condition2)

% condition1='run';condition2='notrun';
% condition1='trial_notrun_highvis';condition2='trial_notrun_lowvis';
% condition1='trial_lowpup_highvis';condition2='trial_lowpup_lowvis';
% condition1='trial_highvis';condition2='trial_lowvis';


end
function findIdentical(condition1, condition2)
trial_notrun_correct_cat_eigenvector=[];trial_notrun_incorrect_cat_eigenvector=[];
animals={'xw','xx','xz','xs','xt','xu'};

for i=1:length(animals)
    animal=char(animals(i)); 
    trial_notrun_correct=load(strcat('X:\Lav\ProcessingDirectory\',animal,'\','network_analysis_corr',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
    trial_notrun_incorrect=load(strcat('X:\Lav\ProcessingDirectory\',animal,'\','network_analysis_corr',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
    trial_notrun_correct_cat_eigenvector=cat(2, trial_notrun_correct_cat_eigenvector,trial_notrun_correct.cent_corr.eigenvector); 
      trial_notrun_incorrect_cat_eigenvector=cat(2, trial_notrun_incorrect_cat_eigenvector,trial_notrun_incorrect.cent_corr.eigenvector); 

end
if any(all(trial_notrun_correct_cat_eigenvector-trial_notrun_incorrect_cat_eigenvector == 0))
disp('Problem with ' )
disp([condition1 ' minus ' condition2]); 
    disp(animals(all(trial_notrun_correct_cat_eigenvector-trial_notrun_incorrect_cat_eigenvector == 0)))
else
    disp('all good'); 
end
end