;  This script will plot the variables output by mcor2.pro
start=-1999
finish=1999



cgplot,lag(start+2000:finish+2000),cc(start+2000:finish+2000,*,1),ticklen=1,xgridstyle=1,ygridstyle=1
cgplot,lag,cc1(*,*,1),/overplot,color='red'
;cgplot,lag,cc2(*,*,1),/overplot,color='grey'
cgplot,lag,cc3(*,*,1),/overplot,color='green'
cgplot,lag,cc4(*,*,1),/overplot,color='blue'

end   
