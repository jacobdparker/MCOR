N=10000


  confidence=.9999
  freedom = 5.

restore, 'moses_super.sav'

;alpha determines the confidence interval of interest
alpha = (1-confidence)^(1/freedom)/2

; these will determine the start and end of the rows used for
; correlation and randomdata generation

s=0
e = 2047


r = randomdata(super_zero,super_minus,super_plus,N,s,e,/long,p_rand=r_p,m_rand=r_m)


r = r(s:e,*)
r_p = r_p(s:e,*)
r_m = r_m(s:e,*)



r_sz = size(r)

;eliminate negative pixel values, this seems to be needed before the
;median can be divided out.

r += make_array(r_sz(1),value=1)#min(r,dimension=1)
r_p += make_array(r_sz(1),value=1)#min(r_p,dimension=1)
r_m += make_array(r_sz(1),value=1)#min(r_m,dimension=1)

;normalize by median before subtraction
r /= make_array(r_sz(1),value=1)#median(r,dimension=1)
r_p /= make_array(r_sz(1),value=1)#median(r_p,dimension=1)
r_m /= make_array(r_sz(1),value=1)#median(r_m,dimension=1)

r_pz = r_p-r
r_mz = r_m-r
;r_pz = r_p


r_mz -= make_array(r_sz(1),value=1)#mean(r_mz,dimension=1)
r_pz -= make_array(r_sz(1),value=1)#mean(r_pz,dimension=1)



lag = [-2047:2047]


;; cc = fltarr(n_elements(lag),r_sz(2))
;; for i = 0,r_sz(2)-1 do begin
;; cc(0,i) = c_correlate(r(*,i),r_pz(*,i),lag)
;; end


r_pad = [r_mz,fltarr(r_sz(1),r_sz(2))]
r_pz_pad = [r_pz,fltarr(r_sz(1),r_sz(2))]


f_r = FFT(r_pad,dimension=1)
f_r_pz = FFT(r_pz_pad,dimension=1)
r_cor = r_sz(1)*2*FFT(f_r*conj(f_r_pz),dimension=1,/inverse)/((fltarr(n_elements(r_pad(*,0)))+1)#sqrt(total(r^2,1)*total(r_pz^2,1)))
r_cor = reverse(real_part(r_cor))
r_cor = shift(r_cor,[-r_sz(1),0])


;time to sort and plot some confidence integrals
for i=0,n_elements(R_cor(*,0))-1 do begin
   R_cor(i,*) = R_cor(i,sort(R_cor(i,*)))
end

;find confdence intervals
k_alpha_r = floor(n_elements(r_cor(0,*))*alpha)
k_1minusalpha_r = floor(n_elements(r_cor(0,*))*(1-alpha))


lag = findgen(n_elements(median(r_cor,dimension=2)))-n_elements(median(r_cor,dimension=2))/2




s_zero = super_zero(s:e,*)
s_plus = super_plus(s:e,*)
s_minus = super_minus(s:e,*)



s_zero /= median(s_zero)
s_plus /= median(s_plus)
s_minus /= median(s_minus)

s_sz = size(s_zero)

s_pz = s_plus-s_zero
s_mz = s_minus-s_zero

s_pz -= make_array(r_sz(1),value=1)#mean(s_pz,dimension=1)
s_mz -= make_array(r_sz(1),value=1)#mean(s_mz,dimension=1)

pz_pad = [s_pz,fltarr(s_sz(1),s_sz(2))]
mz_pad = [s_mz,fltarr(s_sz(1),s_sz(2))]


f_pz = FFT(pz_pad,dimension=1)
f_mz = FFT(mz_pad,dimension=1)
m_cor = r_sz(1)*2*FFT(f_pz*conj(f_mz),dimension=1,/inverse)/((fltarr(n_elements(pz_pad(*,0)))+1)#sqrt(total(s_pz^2,1)*total(s_mz^2,1)))
m_cor = reverse(real_part(m_cor))
m_cor = shift(m_cor,[-s_sz(1),0])

;wishbone_cor = mean(m_cor(*,800:1000),dimension=2)
wishbone_cor =mean(m_cor,dimension=2)

cgwindow,'cgplot',lag,median(r_cor,dimension=2),yr=[-.1,.5],xr=[-1500,1500],xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,$
         xtitle='Lag (pixels)',ytitle='Correlation',wmulti=[0,1,2]
cgwindow,'cgplot',lag,R_cor(*,k_alpha_r),/overplot,color='red',/addcmd
cgwindow,'cgplot',lag,R_cor(*,k_1minusalpha_r),/overplot,color='red',/addcmd

cgwindow,'cgplot',lag,wishbone_cor,/overplot,color='blue',/addcmd

cgwindow,'cgplot',lag,median(r_cor,dimension=2),yr=[0,.5],xr=[-100,100],xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,$
         xtitle='Lag (pixels)',ytitle='Correlation',/addcmd
cgwindow,'cgplot',lag,R_cor(*,k_alpha_r),/overplot,color='red',/addcmd
cgwindow,'cgplot',lag,R_cor(*,k_1minusalpha_r),/overplot,color='red',/addcmd

cgwindow,'cgplot',lag,wishbone_cor,/overplot,color='blue',/addcmd


cgwindow,'cglegend',colors=['black','red','blue'],titles=['Median Correlation',string((confidence)*100,form='(f5.2)')+'% Confidence Interval','Moses Cross Correlation'],Location=[0.15, 0.85],charsize=1,/box,/background,/addcmd


end


  
