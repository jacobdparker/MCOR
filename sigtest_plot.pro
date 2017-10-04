
;also just to be safe run the following
; .r moses_super_exposure
; .r sigtest


;zero_rand is the output random data of sigtest
;super_zero is an output exposure from moses_super_exposure

pro sigtest_plot,zero_rand,super_zero,super_plus,super_minus,alpha,start,auto=auto,power=power,stats=stats,window=window,unitvar=unitvar,disvar=disvar,subtract=subtract

s=150
e = 2042

;Normalize images before substraction
s_zero = super_zero(s:e,*)
s_plus = super_plus(s:e,*)
s_minus = super_minus(s:e,*)

s_zero /= median(s_zero)
s_plus /= median(s_plus)
s_minus /= median(s_minus)

pz = s_plus-s_zero
pz -= (mean(pz,dimension=2)#(fltarr(1024)+1))

z = s_zero
z = z-(mean(z,dimension=2)#(fltarr(1024)+1))
  
if keyword_set(subtract) then z = pz




r_data = zero_rand(start:start+1023,*)



;mean subtract every array   
r_data = r_data -(fltarr(1024)+1)#mean(r_data,dimension=1)
z = super_zero(s:*,*)-(mean(super_zero(s:*,*),dimension=2)#(fltarr(1024)+1))

r_sz = size(r_data)
m_sz = size(z)





if keyword_set(window) then begin
   r_data = (hanning(r_sz(1),alpha=.5)#make_array(r_sz(2),value=1))*r_data
   z = (make_array(m_sz(1),value=1)#hanning(m_sz(2),alpha=.5))*z
end

if keyword_set(disvar) then begin

   var_rand = variance(r_data,dimension=1)
   var_m = variance(z,dimension=2)
   r_data = r_data/(make_array(r_sz(1),value=1)#sqrt(var_rand))
   var_index = round(randomu(seed,r_sz(2))*n_elements(var_m))
   r_data = r_data*(make_array(r_sz(2),value=1)#sqrt(var_m(var_index)))
   

end


   
if keyword_set(unitvar) then begin
   var_rand = variance(r_data,dimension=1)
   var_m = variance(z,dimension=2)
   r_data = r_data/(make_array(r_sz(1),value=1)#sqrt(var_rand))
   z = z/(sqrt(var_m)#make_array(m_sz(2),value=1))
end


if keyword_set(stats) then begin

   m_var = variance(z,dimension=2)
   r_var = variance(r_data,dimension=1)

   hist_m = total(histogram(m_var),/cumulative)/n_elements(m_var)
   hist_r = total(histogram(r_var),/cumulative)/n_elements(r_var)

   window,3
   cgplot,hist_m,color='red',yr=[.8,1]
   cgplot,hist_r,/overplot

end



if keyword_set(auto) then begin
   
;zero pad
r_data_pad = [r_data,fltarr(1024,r_sz(2))]
zero_pad = [[z],[fltarr(n_elements(z(*,0)),1024)]]


f_zero = FFT(zero_pad,dimension=2)
M_auto = 2048*FFT(f_zero*conj(f_zero),dimension=2,/inverse)/(total(z^2,2)#(fltarr(2048)+1))
f_rand = FFT(r_data_pad,dimension=1)
R_auto = 2048*FFT(f_rand*conj(f_rand),dimension=1,/inverse)/((fltarr(2048)+1)#total(r_data^2,1))
M_auto = real_part(m_auto)
r_auto = real_part(r_auto)


;sort for easier statistics
for i=0,n_elements(M_auto(0,*))-1 do begin
   M_auto(*,i) = M_auto(sort(M_auto(*,i)),i)
end
for i=0,n_elements(R_auto(*,0))-1 do begin
   R_auto(i,*) = R_auto(i,sort(R_auto(i,*)))
end

;find confdence intervals
k_alpha_r = round(n_elements(r_auto(0,*))*alpha)
k_1minusalpha_r = round(n_elements(r_auto(0,*))*(1-alpha))
k_alpha_m = round(n_elements(m_auto(*,0))*alpha)
k_1minusalpha_m = round(n_elements(m_auto(*,0))*(1-alpha))


cgwindow,'cgplot',M_auto(k_alpha_m,*),color='red',xr = [0,1024],xtitle='Lag (pixels)',ytitle='Correlations'$
       ,xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1
cgwindow,'cgplot',M_auto(k_1minusalpha_m,*),color='red',/overplot,/addcmd
cgwindow,'cgplot',median(m_auto,dimension=1),color='red',/overplot,/addcmd

cgwindow,'cgplot',R_auto(*,k_alpha_r),/overplot,/addcmd
cgwindow,'cgplot',R_auto(*,k_1minusalpha_r),/overplot,/addcmd
cgwindow,'cgplot',median(r_auto,dimension=2),/overplot,/addcmd
cgwindow,'cglegend',titles=['Moses Columns','Random Data'],colors=['red','black'],location=[.55,.85],/box,/background,/addcmd

endif




if keyword_set(power) then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; compare the data in the fourier domain ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

r_power = abs(FFT(r_data,dimension=1))
m_power = abs(FFT(z,dimension=2))


;sort for easier statistics
for i=0,n_elements(M_power(0,*))-1 do begin
  m_power(*,i) = m_power(sort(M_power(*,i)),i)
end
for i=0,n_elements(R_power(*,0))-1 do begin
   R_power(i,*) = R_power(i,sort(R_power(i,*)))
end

;find confdence intervals
k_alpha_r = round(n_elements(r_power(0,*))*alpha)
k_1minusalpha_r = round(n_elements(r_power(0,*))*(1-alpha))
k_alpha_m = round(n_elements(m_power(*,0))*alpha)
k_1minusalpha_m = round(n_elements(m_power(*,0))*(1-alpha))

cgwindow,'cgplot',M_power(k_alpha_m,*),color='red',/ylog,xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,xtitle='k',ytitle='Power'
cgwindow,'cgplot',M_power(k_1minusalpha_m,*),color='red',/overplot,/addcmd
cgwindow,'cgplot',median(m_power,dimension=1),color='red',/overplot,/addcmd

cgwindow,'cgplot',R_power(*,k_alpha_r),/overplot,/addcmd
cgwindow,'cgplot',R_power(*,k_1minusalpha_r),/overplot,/addcmd
cgwindow,'cgplot',median(r_power,dimension=2),/overplot,/addcmd
cgwindow,'cglegend',titles=['Moses Columns','Random Data'],colors=['red','black'],location=[.5,.85],/box,/background,/addcmd

endif




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
