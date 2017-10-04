

common cor,eit_ch,mcor
string='mcmc_minimization_long.sav'
 restore,string
restore,'mcor.sav'

   mag = [10,1,1,0,0]
      
     
     e = mcor_forward2(mag,tcor=tcor,eit_pz=eit_pz,/plot)
     lag = [-2048:2047]
     cgwindow,'cgplot',lag,tcor,xr=[-1500,1500],yr=[-.5,.5],xtitle='Lag (pixels)',ytitle = 'Correlation'
     cgwindow,'cgplot',lag,mcor,color='red',/overplot,/addcmd
     cgwindow,'cglegend',titles=['MOSES Correlation','Forward Model'],colors=['red','black'],location=[.55,.3],/box,/background,/addcmd
     
     ma = sqrt(total(mcor^2))
     e = 1-(e/ma)
     e *= 100
     error = strtrim(string(e,format='(f5.2)'),1)
     cgwindow,'cgtext',-1350,.4,'Percent Correct = '+error,/addcmd

     pz = (eit_pz-mean(eit_pz))
     pz /= sqrt(total(total(pz^2,2),1))
     
cgwindow,'cgimage',pz,maxvalue=.0005,minvalue=-.0005,ctindex=3,/keep_aspect_ratio,waspect=.5
    cgwindow,'cgpolygon', [0.35, 0.35, 0.65, 0.65, .35], [0.7, 0.97, 0.97, 0.7, 0.7], /NORMAL, COLOR='green',thick='10',/addcmd
cgwindow,'cgpolygon', [0.65, 0.65, 0.95, 0.95, .65], [0.25, 0.5, 0.5, 0.25, 0.25], /NORMAL, COLOR='green',thick='10',/addcmd
cgwindow,'cgpolygon', [0.05, 0.05, 0.12, 0.12, .05], [0.03, 0.2, 0.2, 0.03, 0.03], /NORMAL, COLOR='green',thick='10',/addcmd

                                
s=eit_ch_spectrum(eit_ch,.00025)

  end
