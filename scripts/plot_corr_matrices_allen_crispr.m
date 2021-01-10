function [N, W_spont]=plot_corr_matrices_allen_crispr(outputfolder, animals,statenames, th, parcels_names, parcels_region_labels)

signame='Allen';
%% plot corr matrices
W_spont=nan(23,23,length(statenames), length(animals));
for ai=1:length(animals)
    for state_i=1:length(statenames)
        Wfile=fullfile(outputfolder,animals{ai},[ statenames{state_i} ,signame '_' num2str(th) '.mat']);
        if ~exist(Wfile, 'file')
            continue;
        end
            load(Wfile,'W_corr');
    
     W_spont(:,:,state_i,ai) = process_sim(W_corr, th);
     
        
    end
end
N=sum(~isnan(W_spont(1,1,1,:)));
W_spont=nanmean(W_spont,4);
statenamesstr=statenames;
for si=1:length(statenames)
    statenamesstr{si}(statenames{si}=='_') = ' ';
    
end


figure;
for state_i=1:length(statenames)
    subplot(2,2, state_i);
    imagesc(W_spont(:,:,state_i),[0 1]);colorbar;
    title(['Spont ' statenamesstr{state_i}]);
end
subplot(2,2, 4);
imagesc(W_spont(:,:,end)-W_spont(:,:,1),[-1 1]/10);colorbar;
title(['Spont' statenames{end} '-' statenames{1} ]);
% suptitle(['N=' num2str(N)]);
%%

% 
% figure;
% 
%    
% 
% for state_i=1:length(statenames)
%     subplot(2,2, state_i);
%     G = graph(W_spont(:,:, state_i), parcels_names);
%    p= plot(G, 'LineWidth',G.Edges.Weight*5);
%      p.NodeCData=parcels_region_labels;
%     title(['Spont ' statenamesstr{state_i}]);
% end
%     subplot(2,2, 4);
% G = graph(abs(W_spont(:,:,3)-W_spont(:,:,1)), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*10);
%      p.NodeCData=parcels_region_labels;
% title(['Spont 3-1' ]);
% for state_i=1:length(statenames)
%     subplot(2,2, state_i+4);
%     G = graph(W_c(:,:, state_i), parcels_names);
%   p= plot(G, 'LineWidth',G.Edges.Weight*5);
%      p.NodeCData=parcels_region_labels;
%      title(['Correct ' statenamesstr{state_i}]);
% end
% subplot(2,2, 8);
% G = graph(abs(W_c(:,:,3)-W_c(:,:,1)), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*10);
%      p.NodeCData=parcels_region_labels;
% title(['Correct 3-1' ]);
% for state_i=1:length(statenames)
%      subplot(2,2, state_i+8);
%     G = graph(W_in(:,:, state_i), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*5);
%     p.NodeCData=parcels_region_labels;
%     title(['Inorrect ' statenamesstr{state_i}]);
% end
%      subplot(2,2, 12);
% G = graph(abs(W_in(:,:,3)-W_in(:,:,1)), parcels_names);
% p= plot(G, 'LineWidth',G.Edges.Weight*10);
%      p.NodeCData=parcels_region_labels;
% title(['Inorrect 3-1' ]);
% for state_i=1:length(statenames)
%      subplot(2,2, state_i+12);
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






