restore, 'moses_super.sav'

; Read in all of the eit image from the MOSES launch

read_eit,'/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.184020',index,data
data = iris_prep_despike(data)
eit_prep,index,i304_h,i304,data=data

read_eit,'/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.190014',index,data
data = iris_prep_despike(data)
eit_prep,index,i171_h,i171,data=data


read_eit,'/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.182651',index,data
data = iris_prep_despike(data)
eit_prep,index,i284_h,i284,data=data

read_eit,'/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.183429',index,data
data = iris_prep_despike(data)
eit_prep,index,i195_h,i195,data=data

;form maps for each wavelength
index2map,i304_h,i304,m304
index2map,i171_h,i171,m171
index2map,i195_h,i195,m195
index2map,i284_h,i284,m284


;rotate all maps to match image 13 of MOSES sequence (rotate to
;18:48:19 UT)
m304 = drot_map(m304,.1458,/keep_limb)
m171 = drot_map(m171,-.1858,/keep_limb)
m195 = drot_map(m195,.2433,/keep_limb)
m284 = drot_map(m284,.3706,/keep_limb)

;save a processed cube for moses_eit_linemap
eit_cube = [[[m171.data]],[[m195.data]],[[m284.data]],[[m304.data]]]
save,eit_cube,filename='eit_cube_prepped.sav'

;rebin maps to match MOSES
moses_resolution = .59 ;arcsec per pix
eit_sz = size(i304)
scale = i304_h.cdelt1/moses_resolution
newdim = floor(scale*eit_sz[1])

m171 = rebin_map(m171,newdim,newdim)
m195 = rebin_map(m195,newdim,newdim)
m284 = rebin_map(m284,newdim,newdim)
m304 = rebin_map(m304,newdim,newdim)

;Co-Align and Modify Maps

;pad super_zero
pad_z = fltarr(newdim,newdim)
pad_z(0,0) = super_zero

;correlate both images
cor = convol_fft(m304.data,pad_z,/correlate)
cor = shift(cor,[-(newdim/2-1),-(newdim/2-1)])

m = array_indices(pad_z,where(cor eq max(cor)))

;m=[587,1763] ;result of above corelation to speed things up

m171.data= shift(m171.data,[-m(0),-m(1)])
m284.data= shift(m284.data,[-m(0),-m(1)])
m195.data= shift(m195.data,[-m(0),-m(1)])
m304.data= shift(m304.data,[-m(0),-m(1)])




x = get_map_xp(m304)
x = shift(x,[-m(0),-m(1)])
y = get_map_yp(m304)
y = shift(y,[-m(0),-m(1)])

x = x[0:2047,0:1023]
y = y[0:2047,0:1023]

save,x,y,filename='mosesI_pointing.sav'




;Display EIT images in the appropriate color scale and save them

eit_colors,171,r,g,b
eit_171 = image(alog10(m171.data[0:2047,0:1023]>.1<200),x,y,rgb_table = [[r],[g],[b]],axis_style = 1,buffer=1,title= 'SOHO EIT 171 '+i171_h.date_obs,xtitle = 'X (arcsec)',ytitle = 'Y (arcsec)',xshowtext=0)
eit_colors,284,r,g,b
eit_284 = image(alog10(m284.data[0:2047,0:1023]>.001<50),x,y,rgb_table = [[r],[g],[b]],axis_style = 1,buffer=1,title= 'SOHO EIT 284 '+i284_h.date_obs,xtitle = 'X (arcsec)',ytitle = 'Y (arcsec)',xshowtext=0)
eit_colors,195,r,g,b
eit_195 = image(alog10(m195.data[0:2047,0:1023]>.001<200),x,y,rgb_table = [[r],[g],[b]],axis_style = 1,buffer=1,title= 'SOHO EIT 195 '+i195_h.date_obs,xtitle = 'X (arcsec)',ytitle = 'Y (arcsec)',xshowtext=0)
eit_colors,304,r,g,b
eit_304 = image(alog10(m304.data[0:2047,0:1023]>.001<200),x,y,rgb_table = [[r],[g],[b]],axis_style = 1,buffer=1,title= 'SOHO EIT 304 '+i304_h.date_obs,xtitle = 'X (arcsec)',ytitle = 'Y (arcsec)')



eit_171.save, "../MCOR_Paper/eit_171.eps"
eit_195.save, "../MCOR_Paper/eit_195.eps"
eit_284.save, "../MCOR_Paper/eit_284.eps"
eit_304.save, "../MCOR_Paper/eit_304.eps"

;add new images to idl save file for future use.

i171 = m171.data
i195 = m195.data
i304 = m304.data
i284 = m284.data

save, i171,i171_h,i195,i195_h,i284,i284_h,i304,i304_h,filename='moses_eit.sav'




end
