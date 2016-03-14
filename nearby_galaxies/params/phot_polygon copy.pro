PRO phot_polygon, image_file, ds9reg,out


; PURPOSE:
;     read polygon regions in an image using ds9-format polygon files.
;
; INPUT:
;     image_file - string, name of 2D image 
;
;     ds9reg - string, name of the ds9 region file containing the polygon 
;        information. This file can be generated with ds9 of course.  
;
;           Or, any text file that contains lines starting with 'polygon('
;           (no space or any other characters in front) and end with ')' 
;        can be used.  Lines not starting with 'polygon(' will be igored.  
;        The format should be:
;           polygon(x1,y1,x2,y2,x3,y3,x4,y4...,xn,yn)
;        for polygons with n points.  The polygon loop should not be 
;        closed, i.e., xn not equal to x1, same for y.  There should
;        be no space at all in the line.  xn and yn can be integers or
;        floats.  They should be IMAGE COORDINATES.  RA/Dec will not be
;        accepted.  The points CAN be outside the image.  Multiple polygons
;        can be masked at the same time in one polygon file.
;     pixsize - pixel size of the image
;        
;
;
;


image=mrdfits(image_file,0,hdr)
image=double(image)

;IF n_elements(pix_size) EQ 0 THEN begin
;read,pixsize, prompt='Inter pixel size of the image:'
;GETROT, hdr, rot, cdelt
;pix_size = abs( cdelt)*3600.
;ENDIF

imsize = size(image,/dimen)


phot=dblarr(10)
count = 0 
openr, ds9, ds9reg, /get_lun
string=''



WHILE NOT EOF(ds9) DO BEGIN   ; main loop
readf, ds9, string
IF strmid(string, 0, 8) NE 'polygon(' THEN goto, skip   
; extract coordinate components
coors = strmid(string,8,strlen(string)-9)
coors = float(strsplit(coors,',',/extract))
npoint = n_elements(coors)/2
ind = indgen(npoint)
x = coors[ind*2]-1
y = coors[ind*2+1]-1
merge_vector, x, x[0]  ; to close the loop
merge_vector, y, y[0]

minx = max([min(x),0])
maxx = min([max(x),imsize[0]-1])
miny = max([min(y),0])
maxy = min([max(y),imsize[1]-1])

subimage = image[minx:maxx, miny:maxy]
subsize = size(subimage,/dimen)
xind = lindgen(subsize[0],subsize[1]) mod subsize[0]
yind = lindgen(subsize[0],subsize[1]) / subsize[0]
x = float(x-minx)
y = float(y-miny)

;xind=177.-minx
;yind=303.-miny


; calculate sum of angles
;theta = fltarr(subsize[0],subsize[1])
theta=0.0
FOR i=0,npoint-1 DO BEGIN
   theta1 = atan(y[i]-yind, x[i]-xind)
   theta2 = atan(y[i+1]-yind,x[i+1]-xind)

   dtheta = theta2 - theta1
   A = where(dtheta GT !pi)
   IF total(A) NE -1 THEN dtheta[A] = dtheta[A] - 2*!pi 
   B = where(dtheta LT -!pi)
   IF total(B) NE -1 THEN dtheta[B] = dtheta[B] + 2*!pi
   
   theta = theta + dtheta
;   print,theta1,theta2,dtheta,theta
ENDFOR

A = where(abs(theta) GT !pi) 
IF total(A) NE -1 THEN phot[count]= total(subimage[A],/nan) 
image[minx:maxx, miny:maxy] = subimage
count+=1
skip:
ENDWHILE  ; end of main loop
free_lun, ds9
;stop
;phot=double(phot)

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;calculate surface density
;;;;;;;;;;;;;;;;;;;;;;;;;;;
size_of_arcsec= 3.7815467d      ;pc in distance of M31
size_of_arcsec_kpc= size_of_arcsec*10^(-3d)
area_arcsec_sq= [1455.34d,1649.38d,1455.34d,1649.38d,1649.38d,1552.36d,1455.34d,1455.34d,1649.38d,1649.38d] ;area of regions in arcswc^2
phot_per_arcsec=phot/area_arcsec_sq
area_pc_sq=area_arcsec_sq*(size_of_arcsec^2d) ;area of regions in pc^2
area_kpc_sq=area_arcsec_sq*(size_of_arcsec_kpc^2d) ;area of regions in kpc^2
phot_per_pc_sq=phot/area_pc_sq
phot_per_kpc_sq=phot/area_kpc_sq
;stop
;;;;;;;;;;;;;
;;making_table
;;;;;;;;;;;;;
pub_id=['Bulge','Region 1','Region 2','Region 3','Region 4','Region 5','Region 6','Region 7','Region 8','Region 9']
RAdeg=[10.64583333333333d,10.376708333333333d,11.345208333333332d,10.155708333333331d,10.324416666666664d,10.914874999999999d,10.898833333333332d,10.224916666666665d,10.589999999999998d,10.249999999999998d]
Decdeg=[41.35027777777778d,40.718833333333336d,41.64808333333333d,41.02483333333333d,41.11938888888889d,41.317527777777784d,41.3875d,40.98302777777778d,41.1215d,40.60563888888889d]
ID=['Bulge','irc1','irc2','irc3','irc4','irac5','irc6','irc7','irc8','isocvf']

   data = {Pub_ID:'', RAdeg:0.0, Decdeg:0.0,ID:'',area_arcsec_sq:0.0d,IRAC4_lum:0.0d,IRAC4_lum_per_arcsec_sq:0.0d,IRAC4_lum_per_pc_sq:0.0d,IRAC4_lum_per_kpc_sq:0.0d}
   datas = replicate(data, n_elements(phot))
   datas.Pub_ID = pub_id
   datas.RAdeg = RAdeg
   datas.Decdeg = Decdeg
   datas.ID = ID
   datas.area_arcsec_sq= area_arcsec_sq
   datas.IRAC4_lum=phot
   datas.IRAC4_lum_per_arcsec_sq=phot_per_arcsec
   datas.IRAC4_lum_per_pc_sq=phot_per_pc_sq
   datas.IRAC4_lum_per_kpc_sq=phot_per_kpc_sq
   
   mwrfits,datas, out+'.fits', /create

;stop

END

