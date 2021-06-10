function [regsig, al_coeff] = regress_all_pixels(PixxTime_bl, PixxTime_uv, MAXLEN)

if nargin == 2
    MAXLEN = min(1e4, size(PixxTime_bl, 2));
end
npix = size(PixxTime_bl,1);
al_coeff = zeros(npix,1);
len = min(size(PixxTime_bl,2), size(PixxTime_uv,2));
PixxTime_bl=PixxTime_bl(:,1:len);
PixxTime_uv=PixxTime_uv(:,1:len);
regsig = zeros(size(PixxTime_bl));

for ipix=1:npix
    if mod(ipix,1e3) == 0
        disp(['Regressing pixel ' num2str(ipix) ' of ' num2str(npix)]);
    end
    if any(PixxTime_bl(ipix,:)>0)
        [regsig(ipix,:), al_coeff(ipix)] = regress_pixel(PixxTime_bl(ipix,:), PixxTime_uv(ipix,:), MAXLEN);
    end
    
end