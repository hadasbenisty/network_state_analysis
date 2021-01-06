function main_make_network_figs_2p_lan
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
T = readtable('X:\Lan\FMB208 conditioning imaging data\animal_log_imaging_lt.xlsx');

animals = T.AnimalID;
poptype = T.TargetPop;
dayslist = T.PsychTest;
% i = find(strcmp(animals, 'zb'));
% i(2) = find(strcmp(animals, 'zy'));
% i(3) = find(strcmp(animals, 'zn'));
% i(4) = find(strcmp(animals, 'zi'));
% i(5) = find(strcmp(animals, 'zc'));
% i(6) = find(strcmp(animals, 'zm'));
% i(7) = find(strcmp(animals, 'zk'));
% i(8) = find(strcmp(animals, 'zu'));
% i(9) = find(strcmp(animals, 'zq'));
% i(10) = find(strcmp(animals, 'ze'));
% i(11) = find(strcmp(animals, 'zg'));
% i(12) = find(strcmp(animals, 'zh'));
% i(13) = find(strcmp(animals, 'zv'));

% animals=animals(setdiff(1:length(animals), i));
% poptype=poptype(setdiff(1:length(poptype), i));
% dayslist=dayslist(setdiff(1:length(dayslist), i));

stateslabels = { 'low_pup_q', 'high_pup_q', 'high_pup_l'};
cent_features = {  'eigenvector' 'degree' 'closeness'  'diffmap',  'betweenness' 'pagerank', 'second_eigval'};%};%'eigenvector'
similarity_name = {'pearson_corr',    };%'corr',,  'fullcorr' 'cov''partial_corr'
doover=false;
typeofnrnsvec = {'predictive' 'responsive'  };
%% Fig 4 - network
% plot_centrality_res_per_day(animals, outputfiggolder, stateslabels);
for sim_i = 1:length(similarity_name)
for typenrnsi = 1:length(typeofnrnsvec)
    outputfiggolder = ['X:\Hadas\Meso-imaging\lan\meso_results\figs2p\network_centrality_' similarity_name{sim_i}];
    mkNewDir(outputfiggolder)
    plot_centrality_res_perday(typeofnrnsvec{typenrnsi}, dayslist, poptype, cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels, doover);
    %     plot_centrality_res(poptype, cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels, doover);
    end
end
end

function [correct_states, incorrect_states, spon_states, ...
    correct_slopes, incorrect_slopes, pvals_trials]  = load_centrality_results_perday_predictive(dayslist, cent_features, signame, outputfolder, animals, statenames, isweigtedstr, th)

spon_states = cell(length(animals),1);
correct_states = cell(length(animals),1);
incorrect_states = cell(length(animals),1);
correct_slopes = cell(length(animals),1);
incorrect_slopes = cell(length(animals),1);
pvals_trials = cell(length(animals),1);
valid = cell(length(animals),1);
for i=1:length(animals)
    animal=animals{i};
    currdayslist = eval(dayslist{i});
    valid{i} = zeros(length(currdayslist), length(statenames));
    for day_i = 1:length(currdayslist)
        for state_i = 1:length(statenames)
            if th==-1
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
            else
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            if exist(corfile, 'file')&&...
                    exist(incorfile, 'file' )&&...
                    exist(spontfile, 'file')
                load(corfile, ...
                    'slopes_corr');
                load(incorfile, ...
                    'slopes_incorr');
                
                if size(slopes_corr,2)==1||size(slopes_incorr,2)==1
                    
                    continue;
                end
                valid{i}(day_i, state_i) = 1;
            end
        end
    end
end
for i=1:length(animals)
    animal=animals{i};
    currdayslist = eval(dayslist{i});
    if isempty(valid{i})
        continue;
    end
    for day_i = 1:length(currdayslist)
        if ~all(valid{i}(day_i,:)==1)
            continue;
        end
        corr_slopes=[];incorr_slopes=[];pvalss=[];
        for state_i = 1:length(statenames)
            % trial correct
            if th==-1
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
            else
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            
            load(corfile, 'slopes_corr');
            corr_slopes(:, state_i) = nanmean(slopes_corr,2);
            load(incorfile,'slopes_incorr');
            incorr_slopes(:, state_i) = nanmean(slopes_incorr,2);
            [~, pvalss(:,state_i)] = ttest2(slopes_corr', slopes_incorr', 'tail', 'right');
        end
        correct_slopes{i} = cat(1, correct_slopes{i}, corr_slopes);
        incorrect_slopes{i} = cat(1, incorrect_slopes{i}, incorr_slopes);
        pvals_trials{i} = cat(1, pvals_trials{i}, pvalss);
        centmat = nan(size(slopes_incorr,1), length(statenames),length(cent_features));
        for state_i = 1:length(statenames)
            if th==-1
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
            else
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            
            load(corfile, ['cent_corr_' isweigtedstr ]);
            centvals = eval(['cent_corr_' isweigtedstr ]);
            
            for cent_i = 1:length(cent_features)
                centmat(:, state_i, cent_i) = centvals.(cent_features{cent_i});
            end
        end
        correct_states{i} = cat(1, correct_states{i}, centmat);
        centmat = nan(size(centvals.eigenvector,1), length(statenames),length(cent_features));
        
        for state_i = 1:length(statenames)
            % trial incorrect
            if th==-1
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
            else
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            
            load(incorfile, ['cent_corr_' isweigtedstr ]);
            
            centvals = eval(['cent_corr_' isweigtedstr ]);
            
            for cent_i = 1:length(cent_features)
                centmat(:, state_i, cent_i)  = centvals.(cent_features{cent_i});
            end
        end
        incorrect_states{i} = cat(1, incorrect_states{i}, centmat);
        centmat = nan(size(centvals.eigenvector,1), length(statenames),length(cent_features));
        
        for state_i = 1:length(statenames)
            if th==-1
                
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
            else
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            load(spontfile, ['cent_corr_' isweigtedstr ] );
            centvals = eval(['cent_corr_' isweigtedstr ]);
            
            
            for cent_i = 1:length(cent_features)
                centmat(:, state_i, cent_i) = centvals.(cent_features{cent_i});
            end
        end
        spon_states{i} = cat(1, spon_states{i}, centmat);
    end
    % making sure that missing data is nan
    
end
end



function [correct_states, incorrect_states, spon_states, ...
    correct_slopes, incorrect_slopes, pvals_trials, correct_slopes_baseline, incorrect_slopes_baseline]  = load_centrality_results_perday_responsive(dayslist, cent_features, signame, outputfolder, animals, statenames, isweigtedstr, th)

spon_states = cell(length(animals),1);
correct_states = cell(length(animals),1);
incorrect_states = cell(length(animals),1);
correct_slopes = cell(length(animals),1);
incorrect_slopes = cell(length(animals),1);
correct_slopes_baseline = cell(length(animals),1);
incorrect_slopes_baseline = cell(length(animals),1);
pvals_trials = cell(length(animals),1);
valid = cell(length(animals),1);
for i=1:length(animals)
    animal=animals{i};
    currdayslist = eval(dayslist{i});
    valid{i} = zeros(length(currdayslist), length(statenames));
    for day_i = 1:length(currdayslist)
        for state_i = 1:length(statenames)
            if th==-1
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
            else
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            if exist(corfile, 'file')&&...
                    exist(incorfile, 'file' )&&...
                    exist(spontfile, 'file')
                valid{i}(day_i, state_i) = 1;
            end
        end
    end
end
for i=1:length(animals)
    animal=animals{i};
    currdayslist = eval(dayslist{i});
    if isempty(valid{i})
        continue;
    end
    for day_i = 1:length(currdayslist)
        if ~all(valid{i}(day_i,:)==1)
            continue;
        end
        corr_slopes=[];incorr_slopes=[];pvalss=[];
        for state_i = 1:length(statenames)
            % trial correct
            if th==-1
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
            else
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            
            load(corfile, 'slopes_corr', 'slopes_corr_baseline');
            corr_slopes(:, state_i) = nanmean(slopes_corr,2);
            load(incorfile,'slopes_incorr', 'slopes_inco_baseline');
            incorr_slopes(:, state_i) = nanmean(slopes_incorr,2);
            
            [~, pvalss(:,state_i)] = ttest([slopes_corr,  slopes_incorr]',[slopes_corr_baseline,  slopes_inco_baseline]' , 'tail','left');
        end
        correct_slopes{i} = cat(1, correct_slopes{i}, corr_slopes);
        incorrect_slopes{i} = cat(1, incorrect_slopes{i}, incorr_slopes);
        pvals_trials{i} = cat(1, pvals_trials{i}, pvalss);
        centmat = nan(size(slopes_incorr,1), length(statenames),length(cent_features));
        for state_i = 1:length(statenames)
            if th==-1
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
            else
                corfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            
            load(corfile, ['cent_corr_' isweigtedstr ]);
            centvals = eval(['cent_corr_' isweigtedstr ]);
            
            for cent_i = 1:length(cent_features)
                centmat(:, state_i, cent_i) = centvals.(cent_features{cent_i});
            end
        end
        correct_states{i} = cat(1, correct_states{i}, centmat);
        centmat = nan(size(centvals.eigenvector,1), length(statenames),length(cent_features));
        
        for state_i = 1:length(statenames)
            % trial incorrect
            if th==-1
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
                
            else
                incorfile = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            
            load(incorfile, ['cent_corr_' isweigtedstr ]);
            
            centvals = eval(['cent_corr_' isweigtedstr ]);
            
            for cent_i = 1:length(cent_features)
                centmat(:, state_i, cent_i)  = centvals.(cent_features{cent_i});
            end
        end
        incorrect_states{i} = cat(1, incorrect_states{i}, centmat);
        centmat = nan(size(centvals.eigenvector,1), length(statenames),length(cent_features));
        
        for state_i = 1:length(statenames)
            if th==-1
                
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_ththird.mat']);
            else
                spontfile=fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_D' num2str(currdayslist(day_i)) '_' num2str(th) '.mat']);
            end
            load(spontfile, ['cent_corr_' isweigtedstr ] );
            centvals = eval(['cent_corr_' isweigtedstr ]);
            
            
            for cent_i = 1:length(cent_features)
                centmat(:, state_i, cent_i) = centvals.(cent_features{cent_i});
            end
        end
        spon_states{i} = cat(1, spon_states{i}, centmat);
    end
    % making sure that missing data is nan
    
end
end

function [correct_states, incorrect_states, spon_states]  = load_centrality_results(cent_features, signame, outputfolder, animals, statenames, isweigtedstr, th)

spon_states = cell(length(animals),1);
correct_states = cell(length(animals),1);
incorrect_states = cell(length(animals),1);

for i=1:length(animals)
    animal=animals{i};
    for state_i = 1:length(statenames)
        if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_' num2str(th) '.mat']), 'file')
            
            % trial correct
            load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_' num2str(th)]), ...
                ['cent_corr_' isweigtedstr ]);
            centvals = eval(['cent_corr_' isweigtedstr ]);
            xx=[];
            for cent_i = 1:length(cent_features)
                xx(:, cent_i) = centvals.(cent_features{cent_i});
                
                correct_states{i}(:, state_i, cent_i) = centvals.(cent_features{cent_i});
                
            end
        end
        if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_' num2str(th) '.mat']), 'file' )
            % trial incorrect
            load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_' num2str(th)]), ...
                ['cent_corr_' isweigtedstr ]);
            centvals = eval(['cent_corr_' isweigtedstr ]);
            xx=[];
            for cent_i = 1:length(cent_features)
                xx(:, cent_i) = centvals.(cent_features{cent_i});
                
                incorrect_states{i}(:, state_i, cent_i) = centvals.(cent_features{cent_i});
                
            end
        end
        if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_' num2str(th) '.mat']), 'file')
            % spont
            load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_' num2str(th)]), ...
                ['cent_corr_' isweigtedstr ] );
            centvals = eval(['cent_corr_' isweigtedstr ]);
            
            xx=[];
            for cent_i = 1:length(cent_features)
                xx(:, cent_i) = centvals.(cent_features{cent_i});
                
                spon_states{i}(:, state_i, cent_i) = centvals.(cent_features{cent_i});
                
            end
        end
        
    end
    
end
end


function plot_centrality_res_perday(typeofnrns, dayslist, poptype, cent_features, simname, animals, outputfiggolder, statenames, doover)

outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality_' simname];

signals_names = { 'cells'};%'LSSC''Allen'
isweigtedstr = { 'weighted'};%'notweighted'
for l = 1:length(cent_features)
    mkNewDir(fullfile(outputfiggolder,cent_features{l}));
end
%  accvals = cell(length(animals),1);
% for ai=1:length(animals)
%     for state_i = 1:length(statenames)
%     fileacc = fullfile('X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\behavior_prediction\',...
%         [animals{ai} '_' statenames{state_i} '.mat']);
%     if exist(fileacc, 'file')
%         load(fileacc,'accuracy_mat');
%         accvals{ai}(:, state_i) = mean(accuracy_mat);
%     end
%     end
% end
thvec = [17 27:10:51 7 3 5  ];
thvec=-1;
for sig_i = 1:length(signals_names)
    for isweigted = 1:length(isweigtedstr)
        for th=thvec
            sumfile = fullfile(outputfolder, ['summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'perdays' 'th' num2str(th) typeofnrns '.mat']);
            
            %             sumfile = fullfile(outputfolder, ['summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'perdays' 'th' num2str(th) '.mat']);
            if  ~doover&&exist(sumfile, 'file')
                load(sumfile);
            else
                switch typeofnrns
                    case 'responsive'
                        [correct_states, incorrect_states, spon_states,correct_slopes, incorrect_slopes, pvals_trials, correct_slopes_baseline, incorrect_slopes_baseline] = load_centrality_results_perday_responsive(dayslist, cent_features, signals_names{sig_i}, outputfolder, animals, statenames, isweigtedstr{isweigted}, th);
                        save(sumfile, 'pvals_trials', 'correct_states', 'incorrect_states', 'spon_states', 'correct_slopes', 'incorrect_slopes', 'pvals_trials', 'correct_slopes_baseline', 'incorrect_slopes_baseline');
                    case 'predictive'
                        [correct_states, incorrect_states, spon_states,correct_slopes, incorrect_slopes, pvals_trials] = load_centrality_results_perday_predictive(dayslist, cent_features, signals_names{sig_i}, outputfolder, animals, statenames, isweigtedstr{isweigted}, th);
                        save(sumfile, 'pvals_trials', 'correct_states', 'incorrect_states', 'spon_states', 'correct_slopes', 'incorrect_slopes', 'pvals_trials');
                end
                
                
            end
            
            
            legstr = {'Low Q', 'High Q', 'Loc'};
            populations = unique(poptype);
            diffslopesbypop = cell(length(populations),1);
            diffvec_trial=cell(length(cent_features),length(populations));
            diffvec_spont = cell(length(cent_features),length(populations));
            pvals_trialsFinalbypop= cell(length(populations),1);
             for ai=1:length(correct_states)
                popi = find(strcmp(poptype{ai},populations));
                
                
                for l=find(~strcmp(cent_features, 'second_eigval'))
                    
                    if size(correct_states{ai},2)<3||size(incorrect_states{ai},2)<3||size(incorrect_states{ai},3)<l%||size(accvals{ai},2) < 3
                        continue;
                    end
                    v1 = correct_states{ai}(:, :, l)-incorrect_states{ai}(:, :, l);
                    
                    
                    %                     plot_correct_incorrect_per_state_per_parcels(M, S, 1:size(M,1), statenames(1:size(M,3)))
                    %                     suptitle(cent_features{l});
                    %                     mysave(gcf, fullfile(outputfiggolder,cent_features{l},['trial_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th) animals{ai} '_' poptype{ai}]));
                    v2 = spon_states{ai}(:,3,l)-spon_states{ai}(:,1,l);
                    nonnans = sum(~isnan(v1),2) & ~isnan(v2);
                    diffvec_trial{l,popi} = cat(1,diffvec_trial{l,popi}, v1(nonnans,:));
                    
                    diffvec_spont{l,popi} = cat(1,diffvec_spont{l,popi}, v2(nonnans));
                    if l==1 % do this once
                        v = correct_slopes{ai}-incorrect_slopes{ai};
                         diffslopesbypop{popi} = cat(1,diffslopesbypop{popi}, v(nonnans,:));
                        pvals_trialsFinalbypop{popi} = cat(1,pvals_trialsFinalbypop{popi}, pvals_trials{ai}(nonnans,:));
                    end
                    %                     graph_overlay_allen_3conditions('', '', spon_states{ai}(:,:,1,l),...
                    %                         spon_states{ai}(:,:,2,l),  spon_states{ai}(:,:,3,l),...
                    %                         '',cent_features{l},['3 states ' cent_features{l} ' Centrality '], '',2, legstr);
                    %                     suptitle([animals{ai} '_' poptype{ai}]);
                    %                     mysave(gcf, fullfile(outputfiggolder,cent_features{l}, ['spont_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th) animals{ai} '_' poptype{ai}]));
                    
                    
                    
                end
            end
            for l=find(~strcmp(cent_features, 'second_eigval'))
                figure;ll=1;
                for popi = 1:length(populations)
                    for state_i=1:length(statenames)
                        subplot(3,3,ll);
                        jj = pvals_trialsFinalbypop{popi}(:, state_i) < 0.05;
                        N(popi, state_i) = sum(jj)/length(pvals_trialsFinalbypop{popi}(:, state_i));
                        %                         jj=(diffslopesbypop{popi}(:,state_i))>0.05;
                        %                         plotEmbeddingWithColors([diffvec_spont{l,popi}(jj,:) ...
                        %                             diffvec_trial{l,popi}(jj,state_i)], pvals_trialsFinalbypop{popi}(jj,state_i));
                        %                         plot(diffvec_spont{l,popi}, diffvec_trial{l,popi}(:,state_i),'.');
                        lmodel = fitlm(diffvec_spont{l,popi}(jj,:), diffvec_trial{l,popi}(jj,state_i));
                        plot(lmodel);
                        xlabel('3-1 spont');ylabel('corr-inc');
                        title(sprintf('%s %s R2=%2.2f',populations{popi}, legstr{state_i},  lmodel.Rsquared.Ordinary));
                        leg = get(gca,'Legend');leg.Visible='off';
                        ll=ll+1;
                    end
                end
                suptitle([typeofnrns ' cells']);
                mysave(gcf, fullfile(outputfiggolder,cent_features{l}, ['diffs_',cent_features{l},'_'  isweigtedstr{isweigted}  '_scatter_' signals_names{sig_i} '_th' num2str(th)   'perday_' typeofnrns]));
                figure;bar(N);
                legend(legstr);set(gca,'XTickLabel', populations)
                ylabel(['Fraction of ' typeofnrns ' cells'])
                xlabel('population')
                mysave(gcf, fullfile(outputfiggolder, ['Fraction' typeofnrns]));
                
                figure;ll=1;
                for popi = 1:length(populations)
                    for state_i=1:length(statenames)
                        subplot(3,3,ll);
                        %                         plotEmbeddingWithColors([diffvec_spont{l,popi} ...
                        %                             diffvec_trial{l,popi}(:,state_i)], pvals_trialsFinalbypop{popi}(:,state_i));
                        %              %                         plot(diffvec_spont{l,popi}, diffvec_trial{l,popi}(:,state_i),'.');
                        lmodel = fitlm(diffvec_spont{l,popi}, diffvec_trial{l,popi}(:,state_i));
                        plot(lmodel);
                        leg = get(gca,'Legend');leg.Visible='off';
                        xlabel('3-1 spont');ylabel('corr-inc');title(sprintf('%s %s R2=%2.2f',populations{popi}, legstr{state_i},  lmodel.Rsquared.Ordinary));
                        ll=ll+1;
                    end
                end
                suptitle('All cells');
                mysave(gcf, fullfile(outputfiggolder,cent_features{l}, ['diffs_',cent_features{l},'_'  isweigtedstr{isweigted}  '_scatter_' signals_names{sig_i} '_th' num2str(th)  '_'  'perday' typeofnrns 'All']));
                
                
            end
            
            secondegvali = find(strcmp(cent_features, 'second_eigval'));
            if ~isempty(secondegvali)
                for ai=1:length(animals)
                    if isempty(spon_states{ai}) || size(spon_states{ai},3)<secondegvali
                        vs(:, ai) = nan;
                    else
                        vs(:, ai) = squeeze(spon_states{ai}(1,:,secondegvali));
                    end
                    if isempty(correct_states{ai})|| size(correct_states{ai},3)<secondegvali
                        vc(:, ai) = nan;
                    else
                        vc(:, ai) = correct_states{ai}(1,:,secondegvali);
                    end
                    if isempty(incorrect_states{ai})|| size(incorrect_states{ai},3)<secondegvali
                        vi(:, ai) = nan;
                    else
                        vi(:, ai) = incorrect_states{ai}(1,:,secondegvali);
                    end
                end
                for pi=1:length(populations)
                    ii = strcmp(poptype, populations{pi}) & ~isnan(sum(vs))'...
                        & ~isnan(sum(vc))'& ~isnan(sum(vi))';
                    M = mean(vs(:,ii),2);
                    S=std(vs(:,ii),[],2)/sqrt(sum(ii)-1);
                    figure;subplot(3,1,1);plot_3_bars(M,S,statenames);
                    ylim([0 1]);
                    M(:,1) = mean(vc(:,ii),2);
                    S(:,1)=std(vc(:,ii),[],2)/sqrt(sum(ii)-1);
                    M(:,2) = mean(vi(:,ii),2);
                    S(:,2) = std(vi(:,ii),[],2)/sqrt(sum(ii)-1);
                    title('Spont');
                    subplot(3,1,2);
                    plot_3_bars(M(:,1),S(:,1),statenames);title('correct');ylim([0 1]);
                    subplot(3,1,3);plot_3_bars(M(:,2),S(:,2),statenames);title('incorrect');
                    suptitle(populations{pi});ylim([0 1]);
                    mysave(gcf, fullfile(outputfiggolder,['second_eigval_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)  populations{pi} 'perday' typeofnrns]));
                end
                close all;
            end
            
            
            
        end
    end
end

end

function plot_centrality_res(poptype, cent_features, simname, animals, outputfiggolder, statenames, doover)

outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality_' simname];

signals_names = { 'cells'};%'LSSC''Allen'
isweigtedstr = { 'weighted'};%'notweighted'
for l = 1:length(cent_features)
    mkNewDir(fullfile(outputfiggolder,cent_features{l}));
end
accvals = cell(length(animals),1);
for ai=1:length(animals)
    for state_i = 1:length(statenames)
        fileacc = fullfile('X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\behavior_prediction\',...
            [animals{ai} '_' statenames{state_i} '.mat']);
        if exist(fileacc, 'file')
            load(fileacc,'accuracy_mat');
            accvals{ai}(:, state_i) = mean(accuracy_mat);
        end
    end
end
for sig_i = 1:length(signals_names)
    for isweigted = 1:length(isweigtedstr)
        for th=[ 3 5 7:10:51]
            sumfile = fullfile(outputfolder, ['summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'th' num2str(th) '.mat']);
            if  ~doover&&exist(sumfile, 'file')
                load(sumfile);
            else
                [correct_states, incorrect_states, spon_states] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals, statenames, isweigtedstr{isweigted}, th);
                save(sumfile, 'correct_states', 'incorrect_states', 'spon_states');
                
            end
            
            
            legstr = {'Low Q', 'High Q', 'Loc'};
            populations = unique(poptype);
            accbypop = cell(length(populations),1);
            diffvec_trial=cell(length(cent_features),length(populations));
            diffvec_spont = cell(length(cent_features),length(populations));
            for ai=1:length(correct_states)
                popi = find(strcmp(poptype{ai},populations));
                
                
                for l=1%find(~strcmp(cent_features, 'second_eigval'))
                    M=[];S=[];
                    if size(correct_states{ai},2)<3||size(incorrect_states{ai},2)<3||size(accvals{ai},2) < 3
                        continue;
                    end
                    accbypop{popi} = cat(1, accbypop{popi}, accvals{ai});
                    for k=1:length(statenames)
                        M(:,1,k) = correct_states{ai}(:,k,l);
                        M(:, 2, k) = incorrect_states{ai}(:,k,l);
                        S(:, 1, k) = zeros(size(M(:,1,k)));
                        S(:, 2, k) = zeros(size(M(:,1,k)));
                    end
                    
                    diffvec_trial{l,popi} = cat(1,diffvec_trial{l,popi}, M(:,1,:)-M(:,2,:));
                    plot_correct_incorrect_per_state_per_parcels(M, S, 1:size(M,1), statenames(1:size(M,3)))
                    suptitle(cent_features{l});
                    mysave(gcf, fullfile(outputfiggolder,cent_features{l},['trial_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th) animals{ai} '_' poptype{ai}]));
                    diffvec_spont{l,popi} = cat(1,diffvec_spont{l,popi}, spon_states{ai}(:,3,l)-spon_states{ai}(:,1,l));
                    graph_overlay_allen_3conditions('', '', spon_states{ai}(:,1,l),...
                        spon_states{ai}(:,2,l), spon_states{ai}(:,3,l),...
                        '',cent_features{l},['3 states ' cent_features{l} ' Centrality '], '',2, legstr);
                    suptitle([animals{ai} '_' poptype{ai}]);
                    mysave(gcf, fullfile(outputfiggolder,cent_features{l}, ['spont_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th) animals{ai} '_' poptype{ai}]));
                    
                    
                    
                end
            end
            for l=1%find(~strcmp(cent_features, 'second_eigval'))
                for popi = 1:length(populations)
                    figure;
                    for state_i=1:length(statenames)
                        subplot(2,2,state_i);
                        plotEmbeddingWithColors([diffvec_spont{l,popi} ...
                            squeeze(diffvec_trial{l,popi}(:,:,state_i))], accbypop{popi}(:,state_i));
                        %                 plot(diffvec_spont{l,popi}, squeeze(diffvec_trial{l,popi}(:,:,state_i)),'.');
                        xlabel('3-1 spont');ylabel('corr-inc');title(statenames{state_i});
                    end
                    suptitle(populations{popi});
                    mysave(gcf, fullfile(outputfiggolder,cent_features{l}, ['diffs_',cent_features{l},'_'  isweigtedstr{isweigted}  '_scatter_' signals_names{sig_i} '_th' num2str(th)  '_' populations{popi}]));
                    
                end
            end
            
            secondegvali = find(strcmp(cent_features, 'second_eigval'));
            if ~isempty(secondegvali)
                for ai=1:length(animals)
                    if isempty(spon_states{ai})
                        vs(:, ai) = nan;
                    else
                        vs(:, ai) = spon_states{ai}(1,:,secondegvali);
                    end
                    if isempty(correct_states{ai})
                        vc(:, ai) = nan;
                    else
                        vc(:, ai) = correct_states{ai}(1,:,secondegvali);
                    end
                    if isempty(incorrect_states{ai})
                        vi(:, ai) = nan;
                    else
                        vi(:, ai) = incorrect_states{ai}(1,:,secondegvali);
                    end
                end
                for pi=1:length(populations)
                    ii = strcmp(poptype, populations{pi});
                    M = mean(vs(:,ii),2);
                    S=std(vs(:,ii),[],2)/sqrt(length(ii)-1);
                    figure;subplot(3,1,1);plot_3_bars(M,S,statenames);
                    ylim([0 1]);
                    M(:,1) = mean(vc(:,ii),2);
                    S(:,1)=std(vc(:,ii),[],2)/sqrt(length(ii)-1);
                    M(:,2) = mean(vi(:,ii),2);
                    S(:,2) = std(vi(:,ii),[],2)/sqrt(length(ii)-1);
                    
                    subplot(3,1,2);
                    plot_3_bars(M(:,1),S(:,1),statenames);title('c');ylim([0 1]);
                    subplot(3,1,3);plot_3_bars(M(:,2),S(:,2),statenames);title('i');
                    suptitle(populations{pi});ylim([0 1]);
                    mysave(gcf, fullfile(outputfiggolder,['second_eigval_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)  populations{pi}]));
                end
                close all;
            end
            
            
            disp(1);
        end
    end
end

end


function plot_correct_incorrect_per_3parcels(M, S, parcels_names, statenames,spatialindex)

CondColors = get_3states_colors;
for s=1:length(statenames)
    statenames{s}(statenames{s}=='_') = ' ';
end
figure;
set(gcf,'renderer','painters');
for parcel_i = 1:length(spatialindex)
    subplot(length(spatialindex),1,(parcel_i));
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');hold all;
    h1(1).YData(2:3)=0;
    h1(2).YData(2:3)=0;
    h1(1).FaceColor=CondColors(1,:);
    h1(2).FaceColor=CondColors(1,:);
    h1(2).FaceAlpha=0.4;
    h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');hold all;
    h1(1).YData([1 3])=0;
    h1(2).YData([1 3])=0;
    h1(1).FaceColor=CondColors(2,:);
    h1(2).FaceColor=CondColors(2,:);
    h1(2).FaceAlpha=0.4;h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');
    h1(1).YData([1 2])=0;
    h1(2).YData([1 2])=0;
    h1(1).FaceColor=CondColors(3,:);
    h1(2).FaceColor=CondColors(3,:);
    h1(2).FaceAlpha=0.4;h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    
    title(parcels_names{spatialindex(parcel_i)});
    
    set(gca,'xtick',1:23)
    set(gcf, 'Position',  [1,1, 700,1000]);
    set(gca,'xticklabel',statenames)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);%ylim([0 400]);
end
end


function makepsychbarplot(c50s,animals,statenames,name)
n = length(animals);
mean_mean_across_groups1=nanmean(c50s,1);
std_mean_across_groups1=nanstd(c50s,[],1)./sqrt(n-1);
figure;
CondColors = get_3states_colors;
subplot(1,1,1)
set(gcf,'renderer','Painters')
hold on
for b = 1:3
    bg=bar(b, mean_mean_across_groups1(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
set(gca,'xtick',1:3)
set(gca,'xticklabel',statenames)
h = errorbar(1:3,mean_mean_across_groups1, std_mean_across_groups1,'LineStyle','none','LineWidth',0.5);title(name);
h.Color='k';
set(h, 'marker', 'none');
% mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\','figure_1',name), 'all');
end



