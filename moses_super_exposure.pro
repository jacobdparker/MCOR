;this script attempts to combine all un saturated pixels weighted by
;some equivilent to exposure time (image median?) in order to create
;a single moses image with the greatest signal to noise ratio
;possible in each order.


restore, 'mosesAlignFinal.sav'

cubesz = size(cube_zero)

;indentify good and bad data

fzero = finite(cube_zero)
fplus = finite(cube_plus)
fminus = finite(cube_minus)


;set all nan's to zero so they don't effect the sum
f_cube_zero = cube_zero
f_cube_plus = cube_plus
f_cube_minus = cube_minus

f_cube_zero(where(fzero eq 0))=0
f_cube_plus(where(fplus eq 0))=0
f_cube_minus(where(fminus eq 0))=0


median_zero = fltarr(cubesz(1),cubesz(2),cubesz(3))
median_plus = fltarr(cubesz(1),cubesz(2),cubesz(3))
median_minus = fltarr(cubesz(1),cubesz(2),cubesz(3))



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
badpixtotal = [badpixp,badpixm]

super_plus(badpixtotal)=0
super_zero(badpixtotal)=0
super_minus(badpixtotal)=0

super_plus /= median(super_zero)
super_minus /= median(super_minus)
super_zero /= median(super_zero)

;create subtracted images


pz = super_plus - super_zero
mz = super_minus - super_zero
end
