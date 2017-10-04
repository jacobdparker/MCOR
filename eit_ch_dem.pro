f = file_search('dem/*')
restore,'eit_ch_1_minimized.sav'



ar = read_csv(f(0),types='float')
chole = read_csv(f(1),types='float')
flare = read_csv(f(2),types='float')
prom = read_csv(f(3),types='float')
qs = read_csv(f(4),types='float')
qs_eis = read_csv(f(5),types='float')

dem = [ar.field2*eit_ch.ar.mag,chole.field2*eit_ch.chole.mag,prom.field2*eit_ch.chole.mag,qs.field2*eit_ch.qs.mag,qs_eis.field2*eit_ch.qs_eis.mag]

t = [ar.field1,chole.field1,prom.field1,qs.field1,qs_eis.field1]

h = histogram(t,binsize=.01,min = 4,max = 7,location=lo,reverse_indices=ri)

d = fltarr(n_elements(h))
for i = 0,n_elements(h)-1 do begin
   if ri(i) ne ri(i+1) then begin
      bin_val = ri[ri(i):ri(i+1)-1]
      v= dem(bin_val)
      d(i) = total(v)

    endif
end

j = where(d gt 10)

cgwindow,'cgplot',lo(j),d(j),thick=6,waspect=.6,xtitle='log T (K)',ytitle='log DEM (arbitrary units)',xr=[4.3,6.2],yr=[60,95]

cgwindow,'cgplot',ar.field1,ar.field2/max(ar.field2)*max(d),color='red',psym=-7,/addcmd,/overplot
cgwindow,'cgplot',qs.field1,qs.field2/max(qs.field2)*max(d),/overplot,/addcmd,color='red',psym=-2

;; cgwindow,'cgplot',ar.field1,ar.field2/max(ar.field2)*max(d),color='',/addcmd,/overplot

;; cgwindow,'cgplot',qs.field1,qs.field2/max(qs.field2)*max(d),/overplot,/addcmd,color='red'



cgwindow,'cgplot',flare.field1,flare.field2/max(flare.field2)*max(d),/overplot,/addcmd,color='red'
cgwindow,'cgplot',chole.field1,chole.field2/max(chole.field2)*max(d),/overplot,/addcmd,color='red'
cgwindow,'cgplot',prom.field1,prom.field2/max(prom.field2)*max(d),/overplot,/addcmd,color='red'
cgwindow,'cgplot',qs_eis.field1,qs_eis.field2/max(qs_eis.field2)*max(d),/overplot,/addcmd,color='red'

cgwindow,'cglegend',location=[.4,.9],psym=[0,7,2,0], titles=['Best Fit','Active Region','Quiet Sun','DEM Basis'],colors=['black','red','red','red'],/addcmd,/box,/background



end
