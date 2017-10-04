;before running this script the output of finite_cc_sigtest.pro is
;required through either a run or a restoration of the .sav file

;also just to be safe run the following
; .r moses_super_exposure


;compare distribution of random data to that of real MOSES vertical
;slices.  Going to compare auto correlations of half of the random
;data with auto correlation of vertical slices


r_sz = size(randomdata_zero)

;subtract the mean from the random data and zero pad
;; r_data = randomdata_zero(0:1023,*)
;; r_data = [r_data -(fltarr(1024)+1)#mean(r_data,dimension=1),fltarr(1024,r_sz(2))]
;; zero = super_zero(150:*,*) - (mean(super_zero(150:*,*),dimension=2)#(fltarr(1024)+1))


;; f_zero = FFT([[zero],[fltarr(n_elements(zero(*,0)),1024)]],dimension=2)
;; M_auto = 2048*FFT(f_zero*conj(f_zero),dimension=2,/inverse)/(total(zero^2,2)#(fltarr(2048)+1))
;; f_rand = FFT(r_data,dimension=1)
;; R_auto = 2048*FFT(f_rand*conj(f_rand),dimension=1,/inverse)/((fltarr(2048)+1)#total(r_data^2,1))
;; R_auto = abs(R_auto)
;; M_auto = abs(M_auto)

;; ;sort for easier statistics
;; for i=0,n_elements(M_auto(0,*))-1 do begin
;;    M_auto(*,i) = M_auto(sort(M_auto(*,i)),i)
;; end
;; for i=0,n_elements(R_auto(*,0))-1 do begin
;;    R_auto(i,*) = R_auto(i,sort(R_auto(i,*)))
;; end

;; ;confidence gives the actual interval, k gives the number of data
;; ;points to remove from the ends

;; k = median_confidence(n_elements(R_auto(0,*)),.025,confidence=confidence) 
;; k_m = median_confidence(n_elements(zero(*,0)),.025,confidence=confidence)

;; cgplot,M_auto(n_elements(zero(*,0))/2,0:100),color='red'
;; cgplot,M_auto(n_elements(zero(*,0))/2+k_m,0:100),color='red',/overplot
;; cgplot,M_auto(n_elements(zero(*,0))/2-k_m,0:100),color='red',/overplot

;; cgplot,real_part(R_auto(0:100,r_sz(2)/2)),/overplot
;; cgplot,real_part(R_auto(0:100,r_sz(2)/2+k)),/overplot
;; cgplot,real_part(R_auto(0:100,r_sz(2)/2-k)),/overplot

;; mm=mean(M_auto,dimension=1)
;; rm=mean(R_auto,dimension=2)
;; cgplot,mm(0:200),color='red'
;; cgplot,rm,/overplot

;; cor_sz = size(random_cor)
;; sort_cor = fltarr(cor_sz(1),cor_sz(2))
;; pcor = abs(random_cor)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; compare the data in the fourier domain ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

r_power = FFT(randomdata_zero(0:1023,*),dimension=1)
r_power = r_power*conj(r_power)
r_power = r_power/((fltarr(1024)+1)#r_power(0,*))

m_power = FFT(super_zero(150:*,*),dimension=2)
m_power = m_power*conj(m_power)
m_power = m_power/(m_power(*,0)#(fltarr(1024)+1))

;sort for easier statistics
for i=0,n_elements(M_power(0,*))-1 do begin
  m_power(*,i) = m_power(sort(M_power(*,i)),i)
end
for i=0,n_elements(R_power(*,0))-1 do begin
   R_power(i,*) = R_power(i,sort(R_power(i,*)))
end

k_r = median_confidence(n_elements(R_power(0,*)),.025,confidence=confidence) 
k_m = median_confidence(n_elements(m_power(*,0)),.025,confidence=confidence)

cgplot,M_power(n_elements(m_power(*,0))/2,*),color='red',/ylog
cgplot,M_power(n_elements(m_power(*,0))/2+k_m,*),color='red',/overplot
cgplot,M_power(n_elements(m_power(*,0))/2-k_m,*),color='red',/overplot

cgplot,R_power(*,r_sz(2)/2),/overplot
cgplot,R_power(*,r_sz(2)/2+k),/overplot
cgplot,R_power(*,r_sz(2)/2-k),/overplot






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ploting confidence intervals for correlations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;sort for easier statistics
;; for i=0,cor_sz(1)-1 do begin

;; sort_cor(i,*) = pcor(i,sort(pcor(i,*)))

;; end



;; cgplot,lag,sort_cor(*,n/2),xtitle='Lag',ytitle='Normalized Correlation',title='Median Absolute Correlation at a Given Lag For '+strtrim(string(n,form='(i)'),1)+' Trials'
;; cgplot,lag,sort_cor(*,k),color='red',/overplot
;; cgplot,lag,sort_cor(*,n-k-1),color='red',/overplot
;; cglegend,colors = ['black','red'],titles=['Median Abs. Correlation',string(confidence*100,form='(f5.2)')+'% Confidence Interval'],Location=[0.5, 0.85]


end
