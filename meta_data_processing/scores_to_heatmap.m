function Pvec = scores_to_heatmap(Xmat, tonorm,signame, animal)
switch signame
    case 'Allen'
        Pvec = scores_to_heatmap_allen(Xmat, tonorm);
    case 'LSSC'
        Pvec = scores_to_heatmap_gal(Xmat, animal, tonorm);
    case 'Grid4'
        Pvec = scores_to_heatmap_grid(Xmat, animal, 4);
end

