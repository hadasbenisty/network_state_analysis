function res = get_data_by_state(statesnames, Xa, t_imaging, segments_arousals)
for si = 1:length(statesnames)
    
    eval([statesnames{si} '.Allen = [];']);
    %         eval([statesnames{si} '.Grid4 = [];']);
    eval([statesnames{si} '.t = [];']);
end
for si = 1:length(statesnames)
    if isfield(segments_arousals, statesnames{si})
        data = eval(statesnames{si});
        for seg_i = 1:size(segments_arousals.(statesnames{si}),1)
            
            seg = segments_arousals.(statesnames{si})(seg_i,:);
            ind1 = findClosestDouble(t_imaging, seg(1));
            ind2 = findClosestDouble(t_imaging, seg(2));
            tt=ind1:ind2;
            finalinds = tt(~isnan(sum(Xa(:, tt))));
            data.Allen = cat(2, data.Allen, Xa(:, finalinds));
            %                     data.Grid4 = cat(2, data.Grid4, Xg4(:, finalinds));
            data.t=cat(1,data.t,t_imaging(ind1:ind2));
        end
        
        eval([statesnames{si} '= data;']);
    end
end
for si = 1:length(statesnames)
    res.(statesnames{si}) = eval(statesnames{si});
end