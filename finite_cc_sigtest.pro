;the following must be run inorder to run this program

;; .r moses_super_exposure


;find the power spectrum of vertical columns, average them together,
;interpolate to 2048 wide

pspec_zero = FFT(super_zero(150:*,*),dimension=2)
pspec_zero = pspec_zero*conj(pspec_zero)
median_pspec_zero = median(pspec_zero,dimension=1)
median_pspec_zero = interpol(median_pspec_zero,2048)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; generate white noise power spectrum ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

n=4096*2
randomdata_zero = fltarr(2048,n)


minlag = -2000
maxlag = 2000

lag = findgen(maxlag-minlag)+minlag

random_cor = fltarr(n_elements(lag),n)

for i = 0l,n-1 do begin

rand_index = round(randomu(seed,2048)*2047)
rand_index2 = round(randomu(seed,2048)*2047)
whitepower = fft(randomu(seed,2048))
whitepower2 = fft(randomu(seed,2048))
randomdata_zero(*,i) = fft(whitepower*sqrt(median_pspec_zero),/inverse) 

;random_cor(*,i) = finite_cc(randomdata_zero(*,i),fft(whitepower2*sqrt(median_pspec_zero),/inverse),lag) 

end

;save,lag,randomdata_zero,random_cor,filename='sig_test '+systime()+'.sav'



;test if random data has same power spectrum as mean_pspec

;; av_random_data = total(randomdata_zero,2)/n
;; av_random_pspec = FFT(av_random_data)*conj(FFT(av_random_data))
;; cor_zero = c_correlate (mean_pspec_zero,av_random_pspec,0)





end

     


  
