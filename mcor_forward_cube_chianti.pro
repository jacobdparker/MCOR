min_int = .25                  ;minimun intensity to include in the cube


restore,'moses_eit.sav'
restore,'psfs2.sav'


;grab a section of active region
  
actv_int = median(i304[1200:1450,760:900])

;attempt to normalize intensities to 304

i171 *= actv_int/median(i171[1200:1450,760:900])
i284 *= actv_int/median(i284[1200:1450,760:900])
i195 *= actv_int/median(i195[1200:1450,760:900])

;concatinate eit images

eit = [[[i171]],[[i195]],[[i284]],[[i304]]]

;restore chianti spectrum for active region (as) and quiet sun (qs)

ch_read_fits,'ch_ss_sp_MOSESI_ar.fits',ar
ch_read_fits,'ch_ss_sp_MOSES1_qs.fits',qs

;eit image peak ion_eq temp in order 171 195 284 304 in log T

temps=[6,6.2,6.35,4.7]

;pull relevent info from chianti struct

i = where(ar.lines(*).int gt min_int)
j = where(qs.lines(*).int gt min_int)

ar_elements = ar.lines(i).snote
qs_elements = qs.lines(j).snote


ar_dpix = round(1/.028*(ar.lines(i).wvl-303.7860))  ; caluculate dispersion in the 
qs_dpix = round(1/.028*(qs.lines(j).wvl-303.7860))

ar_int = ar.lines(i).int
qs_int = qs.lines(j).int

ar_wvl = ar.lines(i).wvl
qs_wvl = qs.lines(j).wvl

ar_wvl = ar_wvl(sort(ar_wvl))
qs_wvl = qs_wvl(sort(qs_wvl))


;read in MOSES throughput data
openr,1,'filter1.csv'
filter1 = fltarr(2,300)
readf,1,filter1
close,1

openr,2,'filter2.csv'
filter2 = fltarr(2,300)
readf,2,filter2
close,2

openr,3,'mosesI_grating.csv'
grating = fltarr(5,88)
readf,3,grating
close,3

openr,4,'mosesI_secondary.csv'
secondary = fltarr(4,91)
readf,4,secondary
close,4

;interpolate to wavelength grid from chianti

ar_filter1 = interpol(filter1(1,*),filter1(0,*)*10,ar_wvl)
ar_filter2 = interpol(filter2(1,*),filter2(0,*)*10,ar_wvl)

ar_grating = interpol(grating(1,*),grating(0,*),ar_wvl)

ar_secondary = fltarr(3,n_elements(i))
ar_secondary(0,*) = interpol(secondary(1,*),secondary(0,*),ar_wvl) 
ar_secondary(1,*) = interpol(secondary(2,*),secondary(0,*),ar_wvl)
ar_secondary(2,*) = interpol(secondary(3,*),secondary(0,*),ar_wvl)

qs_filter1 = interpol(filter1(1,*),filter1(0,*)*10,qs_wvl)
qs_filter2 = interpol(filter2(1,*),filter2(0,*)*10,qs_wvl)

qs_grating = interpol(grating(1,*),grating(0,*),qs_wvl)

qs_secondary = fltarr(3,n_elements(j))
qs_secondary(0,*) = interpol(secondary(1,*),secondary(0,*),qs_wvl) 
qs_secondary(1,*) = interpol(secondary(2,*),secondary(0,*),qs_wvl)
qs_secondary(2,*) = interpol(secondary(3,*),secondary(0,*),qs_wvl)


;form combined throughput curves

ar_zero_throughput = 2 * ar_filter1 * ar_filter2 * ar_grating * ar_secondary(0,*)
ar_plus_throughput = ar_filter1 * ar_filter2 * ar_grating * ar_secondary(1,*)
ar_minus_throughput = ar_filter1 * ar_filter2 * ar_grating * ar_secondary(2,*)

qs_zero_throughput = 2 * qs_filter1 * qs_filter2 * qs_grating * qs_secondary(0,*)
qs_plus_throughput = qs_filter1 * qs_filter2 * qs_grating * qs_secondary(1,*)
qs_minus_throughput = qs_filter1 * qs_filter2 * qs_grating * qs_secondary(2,*)

;ar image

pad = fltarr(4096,1024)
sz = size(i)
ar_zero = pad
ar_plus = [pad,pad]
ar_minus = [pad,pad]
for n=0,sz(1)-1 do begin
   k=i(n)
   t = temps-ar.lines(k).tmax
   close_t = where(abs(t) eq min(abs(t)))
   img = eit(*,*,close_t)
   dim = size(img)
   if dim(0) eq 3 then img = mean(img,dimension=3)
   rel_int_zero = ar_zero_throughput(where(ar_wvl eq ar.lines(k).wvl))*ar_int(n)
   rel_int_plus = ar_plus_throughput(where(ar_wvl eq ar.lines(k).wvl))*ar_int(n)
   rel_int_minus = ar_minus_throughput(where(ar_wvl eq ar.lines(k).wvl))*ar_int(n)
   
   ar_zero += img*rel_int_zero(0)
   ar_plus += shift([img*rel_int_plus(0),pad],[-ar_dpix(n),0])
   ar_minus += shift([img*rel_int_minus(0),pad],[ar_dpix(n),0])
   
  
end
ar_zero = ar_zero(0:2047,0:1023)
ar_plus = ar_plus(0:2047,0:1023)
ar_minus = ar_minus(0:2047,0:1023)


;convolve images cubes with MOSES psfs

ar_zero(0,0) = convol_fft(ar_zero(*,*),psfz)
ar_plus(0,0) = convol_fft(ar_plus(*,*),psfp)
ar_minus(0,0) = convol_fft(ar_minus(*,*),psfm)
   

;qs images

sz = size(j)
qs_zero = pad
qs_plus = [pad,pad]
qs_minus = [pad,pad]
for n=0,sz(1)-1 do begin
   k=j(n)
   t = temps-qs.lines(k).tmax
   close_t = where(abs(t) eq min(abs(t)))
   img = eit(*,*,close_t)
   dim = size(img)
   if dim(0) eq 3 then img = mean(img,dimension=3)
   rel_int_zero = qs_zero_throughput(where(qs_wvl eq qs.lines(k).wvl))*qs_int(n)
   rel_int_plus = qs_plus_throughput(where(qs_wvl eq qs.lines(k).wvl))*qs_int(n)
   rel_int_minus = qs_minus_throughput(where(qs_wvl eq qs.lines(k).wvl))*qs_int(n)
   
   qs_zero += img*rel_int_zero(0)
   qs_plus += shift([img*rel_int_plus(0),pad],[-qs_dpix(n),0])
   qs_minus += shift([img*rel_int_minus(0),pad],[qs_dpix(n),0])
   
end
qs_zero = qs_zero(0:2047,0:1023)
qs_plus = qs_plus(0:2047,0:1023)
qs_minus = qs_minus(0:2047,0:1023)


;convolve images cubes with MOSES psfs

qs_zero(0,0) = convol_fft(qs_zero(*,*),psfz)
qs_plus(0,0) = convol_fft(qs_plus(*,*),psfp)
qs_minus(0,0) = convol_fft(qs_minus(*,*),psfm)

eit_ch = {eit_ch,ar_mag:0.0,qs_mag:0.0,ar_zero:ar_zero,ar_plus:ar_plus,ar_minus:ar_minus,qs_zero:qs_zero,qs_plus:qs_plus,qs_minus:qs_minus}
save,eit_ch,filename='eit_ch.sav'



end
