;test the output of finite_cc using two gaussians(one shifted) to help
;interpret real correlation plots











;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; SIMPLE TEST ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x=[-1024:1024]
x=double(x)


zero = exp(-(x^2/(2*(75)^2)))
plus = shift(zero,-2)
minus = shift(zero,2)



cgwindow,'cgplot',x,zero,xticklen=1,xgridstyle=1,xr=[-200,200],wmulti=[0,1,3],xtitle='Lag z',font=.8
cgwindow,'cgplot',x,plus,/overplot,color='red',/addcmd
cgwindow,'cgplot',x,minus,/overplot,color='blue',/addcmd
cgwindow,'cglegend',titles=['Zero Order','Plus Order','Minus Order'],colors=['black','red','blue'],location=[.65,.85],/addcmd

pz = plus-zero
mz = minus-zero

cgwindow,'cgplot',x,pz,xticklen=1,xgridstyle=1,xr=[-200,200],color='red',xtitle='Lag z',font=.8,/addcmd 
cgwindow,'cgplot',x,mz,xticklen=1,xgridstyle=1,color='blue',/overplot,/addcmd
cgwindow,'cglegend',titles=['Plus-Zero Order','Minus-Zero Order'],colors=['red','blue'],location=[.6,.49],/addcmd

sz = size(x)
lag = findgen(3100)-1500
cc = finite_cc(pz,mz,lag)

cgwindow,'cgplot',lag/10,cc,xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,xr=[-200,200],font=.8,xtitle='Lag z',title='Cross Correlation of Plus-Zero and Minus-Zero',/addcmd
;cgwindow,'cglegend',title='Cross Correlation of Plus-Zero and Minus-Zero',location=[.35,.15],/addcmd
end
