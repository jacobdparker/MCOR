;this script attempts to combine all unsaturated pixels weighted by
;some equivalent to exposure time (image median?) in order to create
;a single moses image with the greatest signal to noise ratio
;possible in each order.

;; NOTE:  in order to do this perfectly NANs in saturated regions should be treated as inf and not zero when
;; calculating image median.  Problem is that empty pixels around the edges due to instrument drift are also 
;; marked with NANs.  Treating these pixels as inf would significantly skew the median.  As it stands all NANs
;; are treated as missing data by the IDL median function which produces acceptable results.


restore, 'mosesAlignFinal.sav'

;;; Test that this procedure does the same thing to all three cubes in a couple different ways.

;test = cube_plus
;test1 = cube_minus

;cube_minus = test
;cube_plus = test1

;cube_plus = cube_zero
;cube_minus = cube_zero



; grab the cube size
cubesz = size(cube_zero)

;pick out saturated pixels that are marked as NANs

fzero = finite(cube_zero)
fplus = finite(cube_plus)
fminus = finite(cube_minus)


;form version of cube with zeros in place of saturated pixels
f_cube_zero = cube_zero 
f_cube_plus = cube_plus
f_cube_minus = cube_minus

f_cube_zero(where(fzero eq 0))=0
f_cube_plus(where(fplus eq 0))=0
f_cube_minus(where(fminus eq 0))=0


;create empty arrays the same size as cube_zero
median_zero = fltarr(cubesz(1),cubesz(2),cubesz(3))
median_plus = fltarr(cubesz(1),cubesz(2),cubesz(3))
median_minus = fltarr(cubesz(1),cubesz(2),cubesz(3))



;This for loop simply cycles through each image in the cube and calculates the median, the result is a cube of images where each pixel has the image median as its value
for n=0,cubesz(3)-1 do begin

   median_zero(*,*,n) = median(cube_zero(*,*,n))*fzero(*,*,n)
   median_plus(*,*,n) = median(cube_plus(*,*,n))*fplus(*,*,n)
   median_minus(*,*,n) = median(cube_minus(*,*,n))*fminus(*,*,n)

end


super_zero = total(f_cube_zero,3)/total(median_zero,3)
super_plus = total(f_cube_plus,3)/total(median_plus,3)
super_minus = total(f_cube_minus,3)/total(median_minus,3)

;finally set pixels that apparently never had an information in them
;equal to zero so at least it doesn't appear in a subtracted image.

badpixp=where(finite(super_plus,/nan))
badpixm=where(finite(super_minus,/nan))
badpixtotal = [badpixm,badpixp]

super_plus(badpixtotal)=0
super_zero(badpixtotal)=0
super_minus(badpixtotal)=0


;Finally divide the median from each image so that differences can be taken.
super_plus /= median(super_plus)
super_minus /= median(super_minus)
super_zero /= median(super_zero)


save,super_plus,super_minus,super_zero, filename = 'moses_super.sav'

end
