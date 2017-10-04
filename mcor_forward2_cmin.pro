common cor,eit_ch,mcor

tic
    
start = 'eit_ch_all.sav'
restore,start
restore,'mcor.sav'


he_ii = 1
ar_mag = 0.5
flare_mag = 0.5
qs_mag = 0.5
qs_eis_mag = 0.5
prom_mag =0.5
chole_mag = 0.5

mag = [he_ii,ar_mag,qs_mag,qs_eis_mag,prom_mag,chole_mag,flare_mag]

xbnd = [[0,0,0,0,0,0,0],[10,1,1,1,1,1,1]]
 
gbnd = [[0],[20]]
nobj = 0
gcomp = 'mcor_forward2'

constrained_min,mag,xbnd,gbnd,nobj,gcomp,inform

e = mcor_forward2(mag,tcor=fcor,eit_pz=eit_pz)

save,start,eit_ch,e,mag,fcor,mcor,eit_pz,inform,filename='eit_ch_all_minimized2.sav'
  
toc

end
