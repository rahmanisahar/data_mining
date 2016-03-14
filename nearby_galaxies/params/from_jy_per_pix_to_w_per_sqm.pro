pro from_jy_per_pix_to_w_per_sqm,filename,lambda=lambda,out
;lambda in micron
c= 2.998d8 ;speed of light in vaccume
lambda*=10^(-6d) ; wave in meter
freq=c/lambda

im=mrdfits(filename,0,hd)

im*=10^(-26d) ; from jy to W/m^2/hz

im_f = im * freq ;from W/m^2/hz to W/m^2

sxaddpar,hd,'BUNIT','Wm-2'

writefits,out,im_f,hd

end