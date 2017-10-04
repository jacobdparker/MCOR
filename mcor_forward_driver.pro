common cor,eit_ch,mcor

 restore,'eit_ch.sav'
restore,'mcor.sav'

    
     
     
     e = mcor_forward(mag,tcor=tcor,eit_pz=eit_pz)
     lag = [-2048:2047
     plot,lag,tcor,xr=[-600,600]
     oplot,lag,mcor,color=999

  end
