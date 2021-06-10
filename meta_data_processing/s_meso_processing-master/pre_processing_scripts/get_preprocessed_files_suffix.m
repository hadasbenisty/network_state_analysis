function [fieldname, savesuffix, col] = get_preprocessed_files_suffix(type)
switch type
    case 1
        col=6;% RF mapping
        fieldname = 'RFmap';
        savesuffix='_RFmapping_all';
    case 2
        col=7; % CRF/ORF test before learning days
        fieldname = 'CRF_ORF_before';
        savesuffix='_CRFbefore_all';
    case 3
        col=8; % CS+ CS- Airpuff unpaired
        fieldname = 'Air_CSplus_CS_minus_unpaired';
        savesuffix='_csusunpaired_all';
    case 4
        col=9; % Learning days
        fieldname = 'Layer';
        savesuffix='_learning_all';
    case 5
        col=10; % Learning plateu days
        fieldname = 'Learning_plateu';
        savesuffix='_learningplateau_all';
    case 6
        col=11;  % CSshifted days
        fieldname = 'CS_shift_learning';
        savesuffix='_learningshifted_all';
    case 7
        col=12; % Learning shifted plateu days
        fieldname = 'CS_shift_learning_plateu';
        savesuffix='_learningshiftedplateu_all';
    case 8
        col=13;% Psychtest days
        fieldname = 'Psych_test';
        savesuffix='_psych_all';
    case 9
        col=14;% Psychtest days
        fieldname = 'Psych_test_c50';
        savesuffix='_psychc50_all';
    case 10
        col=15; % CRF/ORF test after learning days
        fieldname = 'CRF_ORF_test_after';
        savesuffix='_CRFafter_all';
    case 11
        col=16; %extinction days
        fieldname = 'extinction';
        savesuffix='_extinction_all';
    case 12
        col=17;
        fieldname = 'Air_CSpl_CSmin_unpaired_after';
        savesuffix='_csusunpaired_after_all';
    case 13
        col=20;
        fieldname = 'Psychtest_full_screen';
        savesuffix='_learning_p1_all';
    case 14
        col=21;
        fieldname = 'Psychtest_Saline_and_drug';
        savesuffix='_learning_p2_all';
    case 15
        col=22;
        fieldname = 'RF_map_Drug';
        savesuffix='_learning_p3_all';
    case 16
        col=21;
        fieldname = 'Psychtest_Saline_and_drug';
        savesuffix='_psychdrug_all';
end