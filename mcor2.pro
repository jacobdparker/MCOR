;; mcor2.pro attempts to do the same thing as mcor.pro which is cross
;; correlate and combine several MOSES images in a way that will
;; identify spectral content of an exposure besides He II. mcor2.pro
;; will do so without ffts.

restore, 'mosesAlignFinal.sav'

;Identify all nans
badpixz=where(finite(cube_zero,/nan))
badpixp=where(finite(cube_plus,/nan))
badpixm=where(finite(cube_minus,/nan))
badpixtotal = [badpixz,badpixp,badpixm]

;Set all bad pix to zero in every image

cube_zero(badpixtotal)=0
cube_plus(badpixtotal)=0
cube_minus(badpixtotal)=0


;Normalize Images by Dividing by the Median

cube_zero = normalize(cube_zero)
cube_minus = normalize(cube_minus)
cube_plus = normalize(cube_plus)

;Generate substracted Image Cubes
mz=cube_minus-cube_zero
pz=cube_plus-cube_zero

;Integrate along the y(non dispersion axis) in each cube

;; cube_zero1D = total(cube_zero,2)
;; cube_minus1D = total(cube_minus,2)
;; cube_plus1D = total(cube_plus,2)

mz1D = total(mz,2) 
pz1D = total(pz,2)

;sectional cross correlations

mz1D1 = total(mz(*,800:1000,*),2) 
pz1D1 = total(pz(*,800:1000,*),2)
mz1D2 = total(mz(*,500:600,*),2) 
pz1D2 = total(pz(*,500:600,*),2)
mz1D3 = total(mz(*,300:500,*),2) 
pz1D3 = total(pz(*,300:500,*),2)
mz1D4 = total(mz(*,0:300,*),2) 
pz1D4 = total(pz(*,0:300,*),2)



;begin generating cross correlation functions for given MOSES image

lag = indgen(4000)-2000
exposure = [0,1,2,3,4,5,6,7,8,9,10]
exposuresize = size(exposure)
lagsize = size(lag)

cc = fltarr(lagsize(1),1,exposuresize(1))

;set up empty sectional cc
ac_p= fltarr(lagsize(1),1,exposuresize(1))
ac_m= fltarr(lagsize(1),1,exposuresize(1))
cc1 = fltarr(lagsize(1),1,exposuresize(1))
cc2 = fltarr(lagsize(1),1,exposuresize(1))
cc3 = fltarr(lagsize(1),1,exposuresize(1))
cc4 = fltarr(lagsize(1),1,exposuresize(1))

for n= 0,exposuresize(1)-1 do begin

cc(*,*,n) = finite_cc(pz1D(*,exposure(n)),mz1D(*,exposure(n)),lag)

;; ;add this section if you want to correlate sections
;; cc1(*,*,n) = finite_cc(pz1D1(*,exposure(n)),mz1D1(*,exposure(n)),lag)
;; cc2(*,*,n) = finite_cc(pz1D2(*,exposure(n)),mz1D2(*,exposure(n)),lag)
;; cc3(*,*,n) = finite_cc(pz1D3(*,exposure(n)),mz1D3(*,exposure(n)),lag)
;; cc4(*,*,n) = finite_cc(pz1D4(*,exposure(n)),mz1D4(*,exposure(n)),lag)

ac_p(*,*,n) = finite_cc(pz1D(*,exposure(n)),pz1D(*,exposure(n)),lag)
ac_m(*,*,n) = finite_cc(mz1D(*,exposure(n)),mz1D(*,exposure(n)),lag)

endfor



end

