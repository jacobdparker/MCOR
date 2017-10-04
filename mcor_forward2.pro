

function mcor_forward2,mag,tcor=tcor,eit_pz=eit_pz,plot=plot
  common cor,eit_ch,mcor
  
  if keyword_set(plot) then mag = 0 else begin

  eit_ch.he_ii.mag = mag(0)
  eit_ch.ar.mag=mag(1)
  eit_ch.qs.mag=mag(2)
  eit_ch.qs_eis.mag = mag(3)
  eit_ch.prom.mag= 0
  eit_ch.chole.mag= 0
  eit_ch.flare.mag = 0

endelse 
  
  
  tz = eit_ch.ar.images(*,*,0)*eit_ch.ar.mag + eit_ch.qs.images(*,*,0)*eit_ch.qs.mag+$
  eit_ch.prom.images(*,*,0)*eit_ch.prom.mag+eit_ch.he_ii.images(*,*,0)*eit_ch.he_ii.mag+$
  eit_ch.chole.images(*,*,0)*eit_ch.chole.mag+eit_ch.flare.images(*,*,0)*eit_ch.flare.mag+eit_ch.qs_eis.images(*,*,0)*eit_ch.qs_eis.mag

  tp = eit_ch.ar.images(*,*,1)*eit_ch.ar.mag + eit_ch.qs.images(*,*,1)*eit_ch.qs.mag+$
  eit_ch.prom.images(*,*,1)*eit_ch.prom.mag+eit_ch.he_ii.images(*,*,1)*eit_ch.he_ii.mag+$
  eit_ch.chole.images(*,*,1)*eit_ch.chole.mag+eit_ch.flare.images(*,*,1)*eit_ch.flare.mag+eit_ch.qs_eis.images(*,*,1)*eit_ch.qs_eis.mag

  tm = eit_ch.ar.images(*,*,2)*eit_ch.ar.mag + eit_ch.qs.images(*,*,2)*eit_ch.qs.mag+$
  eit_ch.prom.images(*,*,2)*eit_ch.prom.mag+eit_ch.he_ii.images(*,*,2)*eit_ch.he_ii.mag+$
  eit_ch.chole.images(*,*,2)*eit_ch.chole.mag+eit_ch.flare.images(*,*,2)*eit_ch.flare.mag+eit_ch.qs_eis.images(*,*,2)*eit_ch.qs_eis.mag
       

  
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
      ;tcor *= max(abs(mcor))/max(abs(tcor))
      
     

      error = sqrt(total((tcor-mcor)^2))
      ;error = max(abs(tcor-mcor))
return, error



end
