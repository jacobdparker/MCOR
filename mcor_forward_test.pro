g =[.0458/2,.0013*2,.0013*2, .0052 ,.0016,.0008,.0032, .09, .0062*2,.0036,.0021*10,.0028,.0045, .0314/4]


restore,'mcor.sav'
  
e = mcor_forward(g,tcor=tcor,eit_pz=eit_pz)

lag = [-2048:2047]


cgplot,lag,tcor,xr=[-100,100]
cgplot,lag,mcor,color='blue',/overplot

end
