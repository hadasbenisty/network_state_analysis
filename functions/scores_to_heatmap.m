function Pmat = scores_to_heatmap(Xvec, tonorm,signame, animal)

switch signame
    case 'LSSC'
        Pmat = scores_to_heatmap_gal(Xvec, animal, tonorm);
    case 'Allen'
         Pmat = scores_to_heatmap_allen(Xvec, tonorm);
    case 'Grid4'
        Pmat = scores_to_heatmap_grid(Xvec, animal, 4);
end
       