function plot_corr_matrices_allen(outputfolder, animals,statenames, th, parcels_names, parcels_region_labels)

signame='Allen';
%% plot corr matrices
W_spont=zeros(23,23,3);
W_c=zeros(23,23,3);
W_in=zeros(23,23,3);
for ai=1:length(animals)
    for state_i=1:length(statenames)
        load(fullfile(outputfolder,[animals{ai} '_',statenames{state_i} ,'_trials_correct_' signame '_W']), 'W_corr_cor');
        load(fullfile(outputfolder,[animals{ai} '_',statenames{state_i} ,'_trials_incorrect_' signame '_W']),'W_corr_inc');
        load(fullfile(outputfolder,[animals{ai} '_',statenames{state_i} ,signame '_W']), 'W_corr');
     W1(:,:,state_i,ai) = process_sim(W_corr_cor, th);
     W2(:,:,state_i,ai) = process_sim(W_corr_inc, th);
     W3(:,:,state_i,ai) = process_sim(W_corr, th);
     
        W_c(:,:, state_i) =  W_c(:,:, state_i)+process_sim(W_corr_cor, th)/length(animals);
        W_in(:,:, state_i) =  W_in(:,:, state_i)+process_sim(W_corr_inc, th)/length(animals);
        W_spont(:,:, state_i) =  W_spont(:,:, state_i)+process_sim(W_corr, th)/length(animals);
        
    end
end
statenamesstr=statenames;
for si=1:length(statenames)
    statenamesstr{si}(statenames{si}=='_') = ' ';
    
end

%%

% 
% figure;
% 
%    
% 
% for state_i=1:length(statenames)
%     subplot(4,4, state_i);
%     G = graph(W_spont(:,:, state_i), parcels_names);
%    p= plot(G, 'LineWidth',G.Edges.Weight*5);
%      p.NodeCData=parcels_region_labels;
%     title(['Spont ' statenamesstr{state_i}]);
% end
%     subplot(4,4, 4);
% G = graph(abs(W_spont(:,:,3)-W_spont(:,:,1)), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*10);
%      p.NodeCData=parcels_region_labels;
% title(['Spont 3-1' ]);
% for state_i=1:length(statenames)
%     subplot(4,4, state_i+4);
%     G = graph(W_c(:,:, state_i), parcels_names);
%   p= plot(G, 'LineWidth',G.Edges.Weight*5);
%      p.NodeCData=parcels_region_labels;
%      title(['Correct ' statenamesstr{state_i}]);
% end
% subplot(4,4, 8);
% G = graph(abs(W_c(:,:,3)-W_c(:,:,1)), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*10);
%      p.NodeCData=parcels_region_labels;
% title(['Correct 3-1' ]);
% for state_i=1:length(statenames)
%      subplot(4,4, state_i+8);
%     G = graph(W_in(:,:, state_i), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*5);
%     p.NodeCData=parcels_region_labels;
%     title(['Inorrect ' statenamesstr{state_i}]);
% end
%      subplot(4,4, 12);
% G = graph(abs(W_in(:,:,3)-W_in(:,:,1)), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*10);
%      p.NodeCData=parcels_region_labels;
% title(['Inorrect 3-1' ]);
% for state_i=1:length(statenames)
%      subplot(4,4, state_i+12);
%      G = graph(abs(W_c(:,:,state_i)-W_in(:,:,state_i)), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*10);
%      p.NodeCData=parcels_region_labels;
%     title([statenamesstr{state_i} ' Corr-Inc']);
% end
% colormap jet
% 
% 
% %%
% figure;
% 
%    
% 
% for state_i=1:length(statenames)
%     subplot(3,3, state_i);
%     G = graph(W_spont(:,:, state_i), parcels_names);
%    p= plot(G, 'LineWidth',G.Edges.Weight*5);
%      p.NodeCData=parcels_region_labels;
%     title(['Spont ' statenamesstr{state_i}]);
%     subplot(3,3, state_i+3);
%     G = graph(W_c(:,:, state_i), parcels_names);
%   p= plot(G, 'LineWidth',G.Edges.Weight*5);
%      p.NodeCData=parcels_region_labels;
%      title(['Correct ' statenamesstr{state_i}]);
%      
%      subplot(3,3, state_i+6);
%     G = graph(W_in(:,:, state_i), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*5);
%     p.NodeCData=parcels_region_labels;
%     title(['Inorrect ' statenamesstr{state_i}]);
%     
% end
% colormap jet






figure;
for state_i=1:length(statenames)
    subplot(4,4, state_i);
    imagesc(W_spont(:,:,state_i),[0 1]);colorbar;
    title(['Spont ' statenamesstr{state_i}]);
end
subplot(4,4, 4);
imagesc(W_spont(:,:,3)-W_spont(:,:,1),[-1 1]/10);colorbar;
title(['Spont 3-1' ]);
for state_i=1:length(statenames)
    subplot(4,4, state_i+4);
    imagesc(W_c(:,:,state_i),[0 1]);colorbar;
    title(['Correct ' statenamesstr{state_i}]);
end

subplot(4,4, 8);
imagesc(W_c(:,:,3)-W_c(:,:,1),[-1 1]/10);colorbar;
title(['Correct 3-1' ]);
for state_i=1:length(statenames)
    subplot(4,4, state_i+8);
    imagesc(W_in(:,:,state_i),[0 1]);colorbar;
    title(['Inorrect ' statenamesstr{state_i}]);
end

subplot(4,4, 12);
imagesc(W_in(:,:,3)-W_in(:,:,1),[-1 1]/10);colorbar;
title(['Inorrect 3-1' ]);
for state_i=1:length(statenames)
    subplot(4,4, state_i+12);
    imagesc(W_c(:,:,state_i)-W_in(:,:,state_i),[-1 1]/30);colorbar;
    title([statenamesstr{state_i} ' Corr-Inc']);
end

