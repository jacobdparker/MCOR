function moses1_eit_linemap,min_int,fits,ions=ions,ar_int=ar_int,ar_wvl=ar_wvl


;inputs
; min_int  is the minimun line intensity to include in the cube
; fits points the the program to the relavant chinati spectrum (a string)
; keyword ions will print out all ions included in the image

  
restore,'moses_eit.sav'
restore,'psfs2.sav'



;restore chianti spectrum 

ch_read_fits,fits,spectra

;read in MOSES throughput data
openr,1,'mosesI_throughput/filter1.csv'
filter1 = fltarr(2,300)
readf,1,filter1
close,1

openr,2,'mosesI_throughput/filter2.csv'
filter2 = fltarr(2,300)
readf,2,filter2
close,2

openr,3,'mosesI_throughput/mosesI_grating.csv'
grating = fltarr(5,88)
readf,3,grating
close,3

openr,4,'mosesI_throughput/mosesI_secondary.csv'
secondary = fltarr(4,91)
readf,4,secondary
close,4

;interpolate to wavelength grid from chianti

;NOTE: ar.lines(*).wvl are not in order.  I hope what the following steps are doing is comparing the uniform transmission files and assigning a through put value for each position in wavelength as listed by Chianti.  Time will tell...

spectra_wvl = spectra.lines(*).wvl

moses_filter1 = interpol(filter1(1,*),filter1(0,*)*10,spectra_wvl)
moses_filter2 = interpol(filter2(1,*),filter2(0,*)*10,spectra_wvl)

moses_grating = interpol(grating(1,*),grating(0,*),spectra_wvl)


moses_secondary = fltarr(3,n_elements(spectra_wvl))
moses_secondary(0,*) = interpol(secondary(1,*),secondary(0,*),spectra_wvl) 
moses_secondary(1,*) = interpol(secondary(2,*),secondary(0,*),spectra_wvl)
moses_secondary(2,*) = interpol(secondary(3,*),secondary(0,*),spectra_wvl)



;form combined throughput curves

moses_zero_throughput = 2 * moses_filter1 * moses_filter2 * moses_grating * moses_secondary(0,*)
moses_plus_throughput = moses_filter1 * moses_filter2 * moses_grating * moses_secondary(1,*)
moses_minus_throughput = moses_filter1 * moses_filter2 * moses_grating * moses_secondary(2,*)


moses_int_zero = spectra.lines.int*moses_zero_throughput
moses_int_plus = spectra.lines.int*moses_plus_throughput
moses_int_minus = spectra.lines.int*moses_minus_throughput

;plot,spectra.lines.wvl,ar_int_zero,psym=2

STOP


;pull relevent info from chianti struct

i = where(ar_int_zero gt min_int)


if total(i) eq -1 then begin

   i = where(ar_int_zero eq max(ar_int_zero))

   print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print,'!! Minimum Intensity Too High.  Only including max line intensity. !!'
   print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

endif

spectra_int = ar_int_zero(i)
spectra_wvl = ar.lines(i).wvl
ions = ar.lines(i).snote

ar_dpix = round(1/.028*(ar.lines(i).wvl-303.7860)) ; caluculate dispersion in the 



;synthetic image

;; pad = fltarr(4096,1024)
;; sz = size(i)
;; ar_zero = pad
;; ar_plus = [pad,pad]
;; ar_minus = [pad,pad]
;; for n=0,sz(1)-1 do begin
;;    k=i(n)
;;    t = [ar.lines(k).tmax,temps]
;;    sort_t = sort(t)
   
   
;;    close = where(sort_t eq 0) -1
;;    closer =  where(sort_t eq 0) +1
;;    if close eq -1 then close=0
;;    if closer eq 5 then closer=4
;;    close = sort_t(close)
;;    closer = sort_t(closer)

;;    alpha = (t(closer)-ar.lines(k).tmax)/(t(closer)-t(close))


   
;;    img = eit(*,*,where(t(closer) eq temps))*(1-alpha(0))+ eit(*,*,where(t(close) eq temps))*alpha(0)
   
;; if closer eq 4 then img = eit(*,*,3)
;; if closer eq 0 then img = eit(*,*,0)   

;;    ar_zero += img*ar_int_zero(k)
;;    ar_plus += shift([img*ar_int_plus(k),pad],[-ar_dpix(n),0])
;;    ar_minus += shift([img*ar_int_minus(k),pad],[ar_dpix(n),0])
   


   
;; end
;; ar_zero = ar_zero(0:2047,0:1023)
;; ar_plus = ar_plus(0:2047,0:1023)
;; ar_minus = ar_minus(0:2047,0:1023)


;convolve images  with MOSES psfs

ar_zero(0,0) = convol_fft(ar_zero(*,*),psfz)
ar_plus(0,0) = convol_fft(ar_plus(*,*),psfp)
ar_minus(0,0) = convol_fft(ar_minus(*,*),psfm)
   

images = [[[ar_zero]],[[ar_plus]],[[ar_minus]]]


return,images

end 
