function contribution_mat = calc_contribution(R2_full, R2_partial, modeled_energy_th)

N = length(R2_full);
contribution_mat = zeros(N);

for n = 1:N
    for k = 1:N
        if n==k || R2_full(n) < modeled_energy_th
            contribution_mat(n, k) = 0;
            continue;
        end
        contribution_mat(n, k) =  max(1-(R2_partial(n, k))/(R2_full(n)),0);
        
%         contribution_mat(n, k) = max( 1-(1-R2_partial(n, k))/(1-R2_full(n)) ,0);
        
    end
end