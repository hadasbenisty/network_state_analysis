function [regsig, al_coeff] = regress_pixel(PixxTime_bl, PixxTime_uv, MAXLEN)
PixxTime_uv=PixxTime_uv(:);
PixxTime_bl=PixxTime_bl(:);
if ~exist('MAXLEN','var')
    MAXLEN = 2e4;
end
tp=min(MAXLEN,length(PixxTime_bl));
al_coeff = regress(PixxTime_bl(1000:tp-1000),PixxTime_uv(1000:tp-1000));%exclude bits of the start and ending when calculating coefficient in case there are weird artifacts there 
regsig=PixxTime_bl-al_coeff*PixxTime_uv;
