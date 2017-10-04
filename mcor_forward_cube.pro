restore,'moses_eit.sav'
restore,'psfs2.sav'


;grab a section of active region
  
actv_int = median(i304[1200:1450,760:900])

;attempt to normalize intensities to 304

i171 *= actv_int/median(i171[1200:1450,760:900])
i284 *= actv_int/median(i284[1200:1450,760:900])
i195 *= actv_int/median(i195[1200:1450,760:900])


;guess from solar spectrum

;g =[.0458,.0013,.0013, .0052 ,.0016,.0008,.0032, .09, .0062,.0036,.0021,.0028,.0045, .0314]

g = fltarr(14)+1


;theoretical spectrum of elements I care about [dpix,relative intensity]
he_ii = [0,1]

fe_xv =[-650,g(0)]
s_xii = [-509,g(1)]
fe_xiv = [-485,g(2)]
ni_xviii = [-391,g(3)]
si_ix = [-365,g(4) ]
si_ix_2 = [-362,g(5)]
si_ix_3 = [-255,g(6)]
si_xi = [-15,g(7)]
mn_xiv = [36,g(8)]
fe_xi = [156,g(9)]
fe_xiii = [278,g(10)]
mg_vii = [504,g(11)]
fe_xiv_2 = [1000,g(12)]
fe_xvi = [1047,g(13)]

spectrum = [[he_ii],[fe_xv],[s_xii],[fe_xiv],[ni_xviii],[si_ix],[si_ix_2],[si_ix_3],[si_xi],[mn_xiv],[fe_xi],[fe_xiii],[mg_vii],[fe_xiv_2],[fe_xvi]]

;assign a representitve eit image to each wavelength
he_ii = i304
fe_xv = i284
s_xii = i284
fe_xiv = i284
ni_xviii = i284
si_ix = i171
si_ix_2 = i171
si_ix_3 = i171
si_xi = i195
mn_xiv = i195
fe_xi = i195
fe_xiii = i284
mg_vii = i171
fe_xiv_2 = i284
fe_xvi = i284

func = [[[he_ii]],[[fe_xv]],[[s_xii]],[[fe_xiv]],[[ni_xviii]],[[si_ix]],[[si_ix_2]],[[si_ix_3]],[[si_xi]],[[mn_xiv]],[[fe_xi]],[[fe_xiii]],[[mg_vii]],[[fe_xiv_2]],[[fe_xvi]]]


tzero=func(*,*,0)
tplus=func(*,*,0)
tminus=func(*,*,0)

i = indgen(n_elements(spectrum(0,*)))
dpix = spectrum(0,i)
rel_int = spectrum(1,i)

sz_s=n_elements(rel_int)

pad = fltarr(2048,1024)
tplus = [tplus,pad]
tminus = [tminus,pad]

;;This Loop will generate a fake moses image cube tzero is the zeroth
;;order,tplus the plus...etc.  

for n=1,sz_s-1 do begin
   zero_scale =func(*,*,n)*rel_int(n)
   tzero=[[[tzero]],[[zero_scale]]]
   plus = shift([zero_scale,pad],[-dpix(n),0])
   minus = shift([zero_scale,pad],[dpix(n),0])
   tminus = [[[tminus]],[[minus]]]
   tplus = [[[tplus]],[[plus]]]
   
end
tzero = tzero(0:2047,0:1023,*)
tplus = tplus(0:2047,0:1023,*)
tminus = tminus(0:2047,0:1023,*)


;convolve images cubes with MOSES psfs

for n=0,sz_s-1 do begin

   tzero(0,0,n) = convol_fft(tzero(*,*,n),psfz)
   tplus(0,0,n) = convol_fft(tplus(*,*,n),psfp)
   tminus(0,0,n) = convol_fft(tminus(*,*,n),psfm)
   
end


;; ;calculate cube median to normalize for subtraction

;; tzero /= median(total(tzero,3))
;; tplus /= median(total(tplus,3))
;; tminus /= median(total(tminus,3))

;; ;compute every possible term in the cross correlation


;; num = n_elements(tzero(0,0,*))
;; cor_cube = fltarr(n_elements(tzero(*,0,0))*2,num,num)
;; ;n_elements(tzero(0,*,0))

;; for j= 0,num-1 do begin
;;    for k = 0,num-1 do begin
;;       s_zero1 = tzero(*,*,j)
;;       s_zero2 = tzero(*,*,k)
;;       s_plus = tplus(*,*,j)
;;       s_minus = tminus(*,*,k)

;;       s_sz = size(s_zero1)

;;       s_pz = s_plus-s_zero1
;;       s_mz = s_minus-s_zero2


      
;;       s_pz -= make_array(s_sz(1),value=1)#mean(s_pz,dimension=1)
;;       s_mz -= make_array(s_sz(1),value=1)#mean(s_mz,dimension=1)

;;       pz_pad = [s_pz,fltarr(s_sz(1),s_sz(2))]
;;       mz_pad = [s_mz,fltarr(s_sz(1),s_sz(2))]


;;       f_pz = FFT(pz_pad,dimension=1)
;;       f_mz = FFT(mz_pad,dimension=1)
;;       test_cor = s_sz(1)*2*FFT(f_pz*conj(f_mz),dimension=1,/inverse)/((fltarr(n_elements(pz_pad(*,0)))+1)#sqrt(total(s_pz^2,1)*total(s_mz^2,1)))
;;       test_cor = reverse(real_part(test_cor))
;;       test_cor = shift(test_cor,[-s_sz(1)+1,0])
;;       test_cor = mean(test_cor,dimension=2)

      
;;       cor_cube(0,j,k) = test_cor
;;    end
;; end




save,tzero,tplus,tminus,filename='eit_cube.sav'

end
