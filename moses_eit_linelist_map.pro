; The goal of this procedure is to turn a dem/pix map generated by
; eit_dem_tool and turn it into a chianti linelist/pix map in order to
; recreate a full disk sun image for each spectral line in the MOSES
; passband. 


;list of EIT images from MOSES I launch
;; f171 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.190014'
;; f195 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.183429'
;; f284 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.182651'
;; f304 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.184020';


;read in EIT images and build a dem map
; 'eit_cube_prepped.sav is produced by moses_eit'
restore, 'eit_cube_prepped.sav'

eit_dem_tool,eit_cube,dem_map,dem_logt


;originally I thought I shouldn't be modeling the low temp, end
;but that is causeing a significant lack of He ii.  Going to include
;it and then hold that part of the DEM Fixed.
;; dem_map = dem_map[*,*,5:-1]
;; dem_logt = dem_logt[5:-1]

;prepare to generate line_lists

wmin = 280
wmax = 340
pressure = 1e15
ioneq_name = '/usr/local/ssw/packages/chianti/dbase/ioneq/chianti.ioneq';obviously matters where you installed chianti

;since ch_synthetic uses EM and not DEM in Isothermal mode we need to
;prepare an array of dt to multiply intensities before
;integrating. I don't know why there is a natural log, but this
;copied from chianti ch_synthetic

temp = 10.^dem_logt
dlnt=ALOG(10.^(dem_logt[1]-dem_logt[0]))
dt=temp * dlnt

;generate an isothermal line list over the range of wavelengths we
;care about, at a constant pressure

ch_synthetic, wmin, wmax, pressure=pressure,logt_isothermal = dem_logt, ioneq_name = ioneq_name,$
              output = line_list, /no_sum_int, all = 0

;need to scale line list to abundance file since
;make_chianti_spec won't do it for /no_sum_int

abund_name = !xuvtop+'/abundance/sun_coronal_1992_feldman.abund' ;another variable internal to chianti
read_abund,abund_name,abund,abund_ref

line_abunds = abund[line_list.lines.iz-1]
line_abunds  = rebin(reform(line_abunds,1,n_elements(line_abunds)),n_elements(line_list.lines[0].int),n_elements(line_abunds))
line_list.lines[*].int = line_list.lines[*].int * line_abunds[*] 




;; ;build line list map
;; dt_len = n_elements(dt)
;; dt = rebin(reform(dt,1,1,dt_len),1024,1024,dt_len)
;; eit_linelist_map = fltarr(1024,1024,n_elements(line_list.lines))

;; for i = 0,n_elements(line_list.lines)-1 do begin
;;    print,i
;;    int_map = rebin(reform(line_list.lines[i].int,1,1,dt_len),1024,1024,dt_len)
;;    eit_linelist_map[0,0,i] = total(int_map*10.^dem_map*dt,3)

;; endfor

meit_ll_map = obj_new('linelistmap')
meit_ll_map->set_dem_map,dem_map
meit_ll_map->set_linelist,line_list

restore, 'mosesI_eit_alignment.sav'  ;apparently this file just store lag numbers for aligning to EIT?

foo = meit_ll_map->rebin_dem_map(newdim,newdim,/overwrite)
foo = meit_ll_map->shift_dem_map(m,/overwrite)

save,meit_ll_map,filename='moses_eit_linelist_map.sav'

obj_destroy,meit_ll_map








     



end
