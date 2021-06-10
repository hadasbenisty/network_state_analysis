function outsig = flipUDsigs(insig,R, C)

for k = 1:size(insig,2)
   x = reshape(insig(:,k), R, C);
   x = flipud(x);
   outsig(:,k) = x(:);    
end