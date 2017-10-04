common cor,tzero,tplus,tminus,mcor

tic

restore,'mcor_forward1.sav'

;g =[.01,.0013,.0013, .0052 ,.0016,.0008,.0032, .05, .0062,.0036,.0021,.0028,.0045, .01]

;second guess

;g =[.1,.0013,.0013, .0052 ,.0016,.0008,.0032, .09, .0062,.0036,.0021,.0028,.0045, .1]

;g = third guess

;g = replicate(.01,14)
g *=2

  restore,'eit_cube.sav'
  restore,'mcor.sav'

  


  xbnd = [[fltarr(n_elements(g))],[fltarr(n_elements(g))+1]]
 

  gbnd = [[0],[2]]
  nobj = 0
  gcomp = 'mcor_forward'

  constrained_min,g,xbnd,gbnd,nobj,gcomp,inform

  e = mcor_forward(g,tcor=fcor,eit_pz=eit_pz)

 save,e,g,fcor,mcor,eit_pz,inform,filename='mcor_forward1_comp.sav'
  
toc

end
