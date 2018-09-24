;; restore,'mosesAlignFinal.sav'

;;   plus = cube_plus[*,*,0]
;;   plus[where(finite(plus,/NAN))]=0
;;   plus /= median(plus)
;;   plus_power = abs(FFT(plus,dimension=2))

;;   zero = cube_zero[*,*,0]
;;   zero[where(finite(zero,/NAN))]=0
;;   zero_power = abs(FFT(zero,dimension=2))

;;   minus = cube_minus[*,*,0]
;;   minus[where(finite(minus,/NAN))]=0
;;   minus_power = abs(FFT(minus,dimension=2))

 

;;   ;sort 
;;   for i=0,n_elements(minus_power(0,*))-1 do minus_power(*,i) = minus_power(sort(minus_power(*,i)),i)
;;   for i=0,n_elements(plus_power(0,*))-1 do plus_power(*,i) = plus_power(sort(plus_power(*,i)),i)
;;   for i=0,n_elements(zero_power(0,*))-1 do zero_power(*,i) = zero_power(sort(zero_power(*,i)),i)


;;   j = 100
;;   p = plot( plus_power[j,*],/ylog)
;;   p = plot( minus_power[j,*],/overplot,color = 'red',/ylog)
;;   p = plot( zero_power[j,*],/overplot,color = 'blue',/ylog)







restore,'moses_super.sav'

plus_power = abs(FFT(super_plus,dimension=2))
minus_power = abs(FFT(super_minus,dimension=2))
;; zero_power = abs(FFT(super_zero,dimension=2))

;sort
for i=0,n_elements(plus_power(0,*))-1 do plus_power(*,i) = plus_power(sort(plus_power(*,i)),i)

for i=0,n_elements(minus_power(0,*))-1 do minus_power(*,i) = minus_power(sort(minus_power(*,i)),i)

;; for i=0,n_elements(zero_power(0,*))-1 do zero_power(*,i) = zero_power(sort(zero_power(*,i)),i)


j = 100
;; p = plot( plus_power[j,*],/ylog)
;; p = plot( minus_power[j,*],/overplot,color = 'red',/ylog)
;; p = plot( zero_power[j,*],/overplot,color = 'blue',/ylog)

p = plot(mean(plus_power[*,*],dimension=1),/ylog)
p = plot(mean(minus_power[*,*],dimension=1),/overplot,color = 'red',/ylog)
p = plot(mean(zero_power[*,*],dimension=1),/overplot,color = 'blue',/ylog)

;fourier filter orders in area of interesting power spectra
p_sz = size(super_plus)
w_plus = (make_array(p_sz(1),value=1)#hanning(p_sz(2),alpha=.5))*super_plus
fft_plus = fft(w_plus,dimension=2)
spectral_window = hanning(180-50,alpha=.5)
big_window = complexarr(p_sz[2])
big_window[50] = spectral_window
big_window[1023-49-(180-50)-1]= spectral_window

real_big_window = rebin(real_part(big_window),p_sz[2],p_sz[1])
real_big_window = rotate(real_big_window,1)
imag_big_window = rebin(imaginary(big_window),p_sz[2],p_sz[1])
imag_big_window = rotate(imag_big_window,1)

big_window = complex(real_big_window,imag_big_window)

filter_plus = real_part(fft(fft_plus*big_window,/inverse,dimension = 2))
;filter_plus = fft(fft_plus*big_window,/inverse,dimension = 2)


w_minus = (make_array(p_sz(1),value=1)#hanning(p_sz(2),alpha=.5))*super_minus
fft_minus = fft(w_minus,dimension=2)

filter_minus = real_part(fft(fft_minus*big_window,/inverse,dimension = 2))




          

end

