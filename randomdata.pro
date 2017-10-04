;the following must be run inorder to run this program

;; .r moses_super_exposure

function randomdata,super_zero,super_minus,super_plus,N,s,e,long=long,p_rand=p_rand,m_rand=m_rand


;Normalize images before substraction
s_zero = super_zero(s:e,*)
s_plus = super_plus(s:e,*)
s_minus = super_minus(s:e,*)

s_zero /= median(s_zero)
s_plus /= median(s_plus)
s_minus /= median(s_minus)

p = s_plus
p -= (mean(p,dimension=2)#(fltarr(1024)+1))

m = s_minus
m -= (mean(m,dimension=2)#(fltarr(1024)+1))

z = s_zero
z = z-(mean(z,dimension=2)#(fltarr(1024)+1))


;attempt to generate random data using a windowed data set
m_sz = size(z)
z = (make_array(m_sz(1),value=1)#hanning(m_sz(2),alpha=.5))*z
p = (make_array(m_sz(1),value=1)#hanning(m_sz(2),alpha=.5))*p
m = (make_array(m_sz(1),value=1)#hanning(m_sz(2),alpha=.5))*m



f_zero = FFT(z,dimension=2)
f_p = FFT(p,dimension=2)
f_m = FFT(m,dimension=2)
sz = size(f_zero)

if keyword_set(long) then begin

;;;;;; this version will create 2048 coulumns by a simalar procedure
zero_rand= fltarr(2048,N)
p_rand = fltarr(2048,N)
m_rand = fltarr(2048,N)

for i=0,N-1 do begin

rand_column = round(randomu(seed,sz(2)+1)*sz(1))
row = indgen(sz(2)+2)/2
row = row(1:n_elements(row)-1)

;; rand_column_p = round(randomu(seed,sz(2)+1)*sz(1))
;; row_p = indgen(sz(2)+2)/2
;; row_p = row_p(1:n_elements(row_p)-1)

;; rand_column_m = round(randomu(seed,sz(2)+1)*sz(1))
;; row_m = indgen(sz(2)+2)/2
;; row_m = row_m(1:n_elements(row_m)-1)

f_rand_m = [f_m(rand_column,row)]
f_rand_p = [f_p(rand_column,row)]
f_rand = [f_zero(rand_column,row)]

;; f_rand_m = [f_m(rand_column_m,row_m)]
;; f_rand_p = [f_p(rand_column_p,row_p)]
;; f_rand = [f_zero(rand_column,row)]

;give a kick to the one Nyquist freqency that shouldn't be one
;; f_rand(1023)=f_rand(1023)*exp(complex(0,1)*randomu(seed)*2*!pi)
;; f_rand_p(1023)=f_rand_p(1023)*exp(complex(0,1)*randomu(seed)*2*!pi)
;; f_rand_m(1023)=f_rand_m(1023)*exp(complex(0,1)*randomu(seed)*2*!pi)

;give a random phase kick to every element
phase = exp(complex(0,1)*randomu(seed,1023)*2*!pi)
f_rand(1)=f_rand(1:1023)*phase

;; phase_p = exp(complex(0,1)*randomu(seed,1023)*2*!pi)
;; f_rand_p(1)=f_rand_p(1:1023)*phase_p

;; phase_m = exp(complex(0,1)*randomu(seed,1023)*2*!pi)
;; f_rand_m(1)=f_rand_m(1:1023)*phase_m

f_rand_p(1)=f_rand_p(1:1023)*phase
f_rand_m(1)=f_rand_m(1:1023)*phase

;reflect about nyquist frequency and chop off the dc term
f_rand = [f_rand,reverse(conj(f_rand(1:n_elements(real_part(f_rand))-2)))]
f_rand_p = [f_rand_p,reverse(conj(f_rand_p(1:n_elements(real_part(f_rand_p))-2)))]
f_rand_m = [f_rand_m,reverse(conj(f_rand_m(1:n_elements(real_part(f_rand_m))-2)))]

zero_rand(0,i) = fft(f_rand,/inverse)
p_rand(0,i) = fft(f_rand_p,/inverse)
m_rand(0,i) = fft(f_rand_p,/inverse)
end

endif else begin

;Test our method by generating 1024 long random coulmns

zero_rand= fltarr(1024,N)
for i=0,N-1 do begin

rand_column = round(randomu(seed,sz(2)/2+1)*sz(1))
row = indgen(sz(2)/2+1)

f_rand = [f_zero(rand_column,row)]

;give a kick to the one Nyquist freqency that shouldn't be one
f_rand(1023)=f_rand(1023)*exp(complex(0,1)*randomu(seed)*2*!pi)

;give a random phase kick to every element
phase = exp(complex(0,1)*randomu(seed,511)*2*!pi)
f_rand(1)=f_rand(1:511)*phase

;reflect about nyquist frequency and chop off the dc term
f_rand = [f_rand,reverse(conj(f_rand(1:n_elements(real_part(f_rand))-2)))]

;invert the randomly generated fourier transform 

zero_rand(0,i) = fft(f_rand,/inverse)
;; zero_rand(0,i) = zero_rand(*,i)-mean(zero_rand(*,i))
;; var_rand = variance(zero_rand(*,i))
;; var_index = round(randomu(seed)*n_elements(z(s:*,0)))
;; zero_rand(0,i) = zero_rand(*,i)/sqrt(var_rand)*sqrt(var(var_index))
end

end


return,zero_rand

end

     


  
