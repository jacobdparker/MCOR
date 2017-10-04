


;test the output of finite_cc using two gaussians(one shifted) to help
;interpret real correlation plots

function mcor_forward,g,tcor=tcor,eit_pz=eit_pz
  common cor,tzero,tplus,tminus,mcor

  
  ;; g1 = g#g
  ;; temp = cor_cube
  ;; for i =1,n_elements(g) do begin
  ;;    for j = 1,n_elements(g) do begin
  ;;       temp(0,i,j) =  g1(i-1,j-1)*temp(*,i,j)
  ;;    end
     
  ;; end
  tg = [2,g]
  tg = reform(tg,1,1,n_elements(tg))
  tg = rebin(tg,n_elements(tzero(*,0,0)),n_elements(tzero(0,*,0)),n_elements(tg))

  tz = total(tzero*tg,3)
  tm = total(tminus*tg,3)
  tp = total(tplus*tg,3)
  
  
;calculate cube median to normalize for subtraction

  tz /= median(tz)
  tp /= median(tp)
  tm /= median(tm)

;compute cc
  
      s_sz = size(tz)

      s_pz = tp-tz
      s_mz = tm-tz

      eit_pz = s_pz
      
      s_pz -= make_array(s_sz(1),value=1)#mean(s_pz,dimension=1)
      s_mz -= make_array(s_sz(1),value=1)#mean(s_mz,dimension=1)

      pz_pad = [s_pz,fltarr(s_sz(1),s_sz(2))]
      mz_pad = [s_mz,fltarr(s_sz(1),s_sz(2))]


      f_pz = FFT(pz_pad,dimension=1)
      f_mz = FFT(mz_pad,dimension=1)
      test_cor = s_sz(1)*2*FFT(f_pz*conj(f_mz),dimension=1,/inverse)/((fltarr(n_elements(pz_pad(*,0)))+1)#sqrt(total(s_pz^2,1)*total(s_mz^2,1)))
      test_cor = reverse(real_part(test_cor))
      test_cor = shift(test_cor,[-s_sz(1)+1,0])
      tcor = mean(test_cor,dimension=2)
      tcor *= max(mcor)/max(tcor)
      
     



;try fitting just to center and not whole function
     
error = sqrt(total((tcor-mcor)^2))


return, error



end
