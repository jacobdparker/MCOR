function eit_ch_spectrum,eit_ch,min_int


wvl = eit_ch.ar.wvl
int = eit_ch.ar.int*eit_ch.ar.mag
ions = eit_ch.ar.ions

wvl = [wvl,eit_ch.he_ii.wvl]
int = [int,eit_ch.he_ii.int*eit_ch.he_ii.mag]
ions = [ions,eit_ch.he_ii.ions]

wvl = [wvl,eit_ch.qs_eis.wvl]
int = [int,eit_ch.qs_eis.int*eit_ch.qs_eis.mag]
ions = [ions,eit_ch.qs_eis.ions]

wvl = [wvl,eit_ch.qs.wvl]
int = [int,eit_ch.qs.int*eit_ch.qs.mag]
ions = [ions,eit_ch.qs.ions]

wvl = [wvl,eit_ch.prom.wvl]
int = [int,eit_ch.prom.mag*eit_ch.prom.int]
ions = [ions,eit_ch.prom.ions]

wvl = [wvl,eit_ch.chole.wvl]
int = [int,eit_ch.chole.mag*eit_ch.chole.int]
ions = [ions,eit_ch.chole.ions]

wvl = [wvl,eit_ch.flare.wvl]
int = [int,eit_ch.flare.int*eit_ch.flare.mag]
ions = [ions,eit_ch.flare.ions]


i = sort(wvl)
wvl = wvl(i)
ions = ions(i)
int = int(i)

h = histogram(wvl,binsize=.1,min=283,max=336,reverse_indices=ri,locations=lo)

w = fltarr(n_elements(h))
line = strarr(n_elements(h))
for i = 0,n_elements(h)-1 do begin
   if ri(i) ne ri(i+1) then begin
      bin_val = ri[ri(i):ri(i+1)-1]
      w(i) = total(int(bin_val))
      el = ions(bin_val)
      line(i) = el(where(int(bin_val) eq max(int(bin_val))))
   endif


end




;read in MOSES throughput data
openr,1,'mosesI_throughput/filter1.csv'
filter1 = fltarr(2,300)
readf,1,filter1
close,1

openr,2,'mosesI_throughput/filter2.csv'
filter2 = fltarr(2,300)
readf,2,filter2
close,2

openr,3,'mosesI_throughput/mosesI_grating.csv'
grating = fltarr(5,88)
readf,3,grating
close,3

openr,4,'mosesI_throughput/mosesI_secondary.csv'
secondary = fltarr(4,91)
readf,4,secondary
close,4

;interpolate to wavelength grid from chianti

ar_wvl = lo

ar_filter1 = interpol(filter1(1,*),filter1(0,*)*10,ar_wvl)
ar_filter2 = interpol(filter2(1,*),filter2(0,*)*10,ar_wvl)

ar_grating = interpol(grating(1,*),grating(0,*),ar_wvl)


ar_secondary = fltarr(3,n_elements(ar_wvl))
ar_secondary(0,*) = interpol(secondary(1,*),secondary(0,*),ar_wvl) 
ar_secondary(1,*) = interpol(secondary(2,*),secondary(0,*),ar_wvl)
ar_secondary(2,*) = interpol(secondary(3,*),secondary(0,*),ar_wvl)



;form combined throughput curves

ar_zero_throughput = 2 * ar_filter1 * ar_filter2 * ar_grating * ar_secondary(0,*)

w *= max(ar_zero_throughput)/max(w)

cgwindow,'cgplot',lo,w,xtitle='Wavelength($\Angstrom$)',ytitle='Transmission',color='red',yr=[0,.02]

cgwindow,'cgplot',lo,ar_zero_throughput,/overplot,/addcmd
cgwindow,'cglegend',titles=['Spectral Contribution (Arbitrary Units)','MOSES Throughput Curve'],colors=['red','black'],/addcmd

n_lines = total(w gt 0)
n_lines=strtrim(string(n_lines,format='(f5.0)'),1)


cgwindow,'cgtext',310,.1,'# included lines = '+n_lines,/addcmd

j = where(w gt min_int)

for n =0,n_elements(j)-1 do cgwindow,'cgtext',lo(j),w(j),line(j),orientation=90,/data,/addcmd  

tint = total(int)
he = total(int(where(ions eq 'He II   ')))/tint
he_ext = eit_ch.he_ii.int*eit_ch.he_ii.mag
he_ext = he_ext/tint
si= total(int(where(ions eq 'Si XI   ')))/tint


print,'Percent Helium'
print,he*100
print,'Percent He unaccounted for by Chianti'
print,he_ext*100
print,'Percent SI XI'
print,si*100
print,'Percent Other'
print,(1-he-si)*100
end

