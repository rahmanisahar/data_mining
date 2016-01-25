pro from_Mjy_per_sr_to_W_per_sqm, filename,lambda=lambda,out
;lambda in micron
c= 2.998d8 ;speed of light in vaccume
lambda*=10^(-6d) ; wave in meter
freq=c/lambda


im=mrdfits(filename,0,hd)

pix_scale, hd, pix_size ; in degrees
pix_size = pix_size*3600. ; in arc seconds
area_pix= pix_size^2d
sr_per_pix= area_pix*2.3504d*10^(-11d)

im2= im * sr_per_pix ;MJy/pix

im_j=im2*1d6 ;Jy/pix

im4= im_j*10^(-26d) ;W/m2/Hz

im_f=im4 * freq ;W/m2

sxaddpar,hd,'BUNIT','W/m2'

writefits,out,im_f,hd

end
