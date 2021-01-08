function plot_bars_by_conditions(spon_states, str, parcels_names, legstr, curtype)

if length(legstr)==3
     num = sum(~isnan(spon_states(1,:,1)));
graph_overlay_allen_3conditions('', '', spon_states(:,:,1),...
                        spon_states(:,:,2), spon_states(:,:,3),...
                        '',str,[curtype  ' '  str ' Centrality '], parcels_names,num, legstr);
else
  graph_overlay_allen_2conditions('', '', spon_states(:,:,1),...
                        spon_states(:,:,2),...
                        '',str,[curtype  ' '  str ' Centrality '], parcels_names, legstr);


end