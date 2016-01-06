
PRO merge_vector, v1, v2

; merge two vectors (1-D) into one
; v1 will be replaced with the new merged vector
; 



n1 = n_elements(v1)
n2 = n_elements(v2)

v3 = make_array(n1+n2, type=size(v1,/type))

v3[0:n1-1] = v1
v3[n1:n1+n2-1] = v2
v1 = v3


END


pro phot_for_mining, file, out

;dirin='~/Desktop/project/data_mining/m31/images/'
dirout='~/Desktop/project/data_mining/m31/fits_tables/'

image=mrdfits(file,0,hdr)
image=double(image)

imsize = size(image,/dimen)

ra1=[10.642846, 10.374729, 11.3433, 10.153179, 10.321883, 10.912417, 10.896379,10.587492, 10.247504, 10.222404]
dec1=[41.357708, 40.726258, 41.655525, 41.032617, 41.127169, 41.325325, 41.395297, 41.129289, 40.613422, 40.990814]
ra2=[10.635521,10.366492,11.334863,10.144958,10.31365,10.904092,10.888038,10.579233,10.239317,10.214175]
dec2=[41.352303,40.720147,41.649478,41.026461,41.121017,41.319222,41.389203,41.123156,40.607286,40.984675]
ra3=[10.648479,10.378825,11.347242,10.157917,10.326621,10.916996,10.900942,10.592162,10.252154,10.227092]
dec3=[41.342408,40.7106,41.639833,41.016617,41.111167,41.309294,41.379264,41.113275,40.597411,40.974806]
ra4=[10.6558,10.387063,11.355679,10.166133,10.334854,10.925321,10.909283,10.600421,10.260346,10.235325]
dec4=[41.347811,40.716711,41.645881,41.022769,41.117319,41.315394,41.385358,41.119408,40.603547,40.980944]

adxy,hdr,ra1,dec1,x1,y1

adxy,hdr,ra2,dec2,x2,y2

adxy,hdr,ra3,dec3,x3,y3

adxy,hdr,ra4,dec4,x4,y4

write_csv, '~/Desktop/project/data_mining/M31/ds9_reg/ds9_xy.reg','ploygon('+STRTRIM(x1, 1) , y1 , x2, y2 , x3, y3, x4, STRTRIM(y4, 1)+')'

;get_lun, unit
openw, unit, '~/Desktop/project/data_mining/M31/ds9_reg/ds9_xy.reg',/Get_LUN, WIDTH=250
printf, unit, '# Region file format: DS9 version 4.1'
printf, unit, 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
printf, unit, 'image'
FOR i=0, n_elements(x1)-1 DO printf, unit, 'polygon(',x1[i],',',y1[i],',',x2[i],',',y2[i],',',x3[i],',',y3[i],',', x4[i],',',y4[i],')'
close, unit
free_lun, unit

ds9reg='~/Desktop/project/data_mining/M31/ds9_reg/ds9_xy.reg'
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
;phot_per_arcsec=phot/area_arcsec_sq
;area_pc_sq=area_arcsec_sq*(size_of_arcsec^2d) ;area of regions in pc^2
;area_kpc_sq=area_arcsec_sq*(size_of_arcsec_kpc^2d) ;area of regions in kpc^2
;phot_per_pc_sq=phot/area_pc_sq
;phot_per_kpc_sq=phot/area_kpc_sq
;stop
;;;;;;;;;;;;;
;;making_table
;;;;;;;;;;;;;
pub_id=['Region 10','Region 1','Region 2','Region 3','Region 4','Region 5','Region 6','Region 7','Region 8','Region 9']
RAdeg=[10.64583333333333d,10.376708333333333d,11.345208333333332d,10.155708333333331d,10.324416666666664d,10.914874999999999d,10.898833333333332d,10.224916666666665d,10.589999999999998d,10.249999999999998d]
Decdeg=[41.35027777777778d,40.718833333333336d,41.64808333333333d,41.02483333333333d,41.11938888888889d,41.317527777777784d,41.3875d,40.98302777777778d,41.1215d,40.60563888888889d]
;ID=['Bulge','irc1','irc2','irc3','irc4','irac5','irc6','irc7','irc8','isocvf']

   data = {Pub_ID:'', RAdeg:0.0, Decdeg:0.0, pacs100:0.0d}
   datas = replicate(data, n_elements(phot))
   datas.Pub_ID = pub_id
   datas.RAdeg = RAdeg
   datas.Decdeg = Decdeg
   ;datas.ID = ID
   ;datas.area_arcsec_sq= area_arcsec_sq
   datas.pacs100=phot
   
   mwrfits,datas, dirout+out+'.fits', /create

;stop

END


