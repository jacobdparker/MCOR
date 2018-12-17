;take a linelist map object and build synthetic MOSES plus, minus, and
;zero orders for a given dem basis

restore, "moses_eit_linelist_map.sav"

ll = meit_ll_map->get_linelist()




;Pair down line list since I don't have all day
percent_max_int= .001 ;with new method I bet I can include every line
lines_small = ll.lines[where(total(ll.lines.int,1) gt percent_max_int/100*max(total(ll.lines.int,1)))]


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

;interpolate to wavelength grid from chianti line list

;; ll_wvl = ll.lines(*).wvl
ll_wvl = lines_small.wvl

filter1 = interpol(filter1(1,*),filter1(0,*)*10,ll_wvl)
filter2 = interpol(filter2(1,*),filter2(0,*)*10,ll_wvl)

grating = interpol(grating(1,*),grating(0,*),ll_wvl)


secondary_orders = fltarr(3,n_elements(ll_wvl))
secondary_orders(0,*) = interpol(secondary(1,*),secondary(0,*),ll_wvl) 
secondary_orders(1,*) = interpol(secondary(2,*),secondary(0,*),ll_wvl)
secondary_orders(2,*) = interpol(secondary(3,*),secondary(0,*),ll_wvl)



;form combined throughput curves, the zero order has two filters and
;is therefore squared.  The index in each order should match the index
;in linelist.lines

zero_throughput = (filter1 * filter2)^2 * grating * secondary_orders(0,*) 
plus_throughput = filter1 * filter2 * grating * secondary_orders(1,*)
minus_throughput = filter1 * filter2 * grating * secondary_orders(2,*)


;;;;;;; Note: this method will not work as implemented as the scaling
;;;;;;; of a DEM basis function cannot be done out side of this
;;;;;;; procedure.  I don't think there is any operation that
;;;;;;; will scale an image the same way scaling the DEM in this code does.

;form DEM basis, for now triangle functions in logt/logdem space
;; logt = ll.logt_isothermal
;; logt_sz = n_elements(logt)

;; basis_elements = 7
;; mag = 15

;; x = findgen(logt_sz)
;; x0 = findgen(basis_elements)/(basis_elements-1)*logt_sz-logt_sz/2.


;; dem_basis = fltarr(logt_sz,basis_elements)
;; for i = 0,basis_elements-1 do begin
;;    dem_basis[0,i] = -abs(2./logt_sz*(x-x0[i])-1)*mag
;; end
;; dem_basis +=30

;make empty array for building sythetic images
dem_map_sz = size(meit_ll_map->get_dem_map())

basis_elements = dem_map_sz[3]

moses_synth_cube = fltarr(dem_map_sz[1],dem_map_sz[2],3,basis_elements)


;since ch_synthetic uses EM and not DEM in Isothermal mode we need to
;prepare an array of dt to multiply intensities before
;integrating. I don't know why there is a natural log, but it is
;copied from chianti ch_synthetic

temp = 10.^ll.logt_isothermal
dlnt = ALOG(10.^(ll.logt_isothermal[1]-ll.logt_isothermal[0]))
dt = temp * dlnt


;-- rebin and reform to multiply against dem map

;; dt_len = n_elements(dt)
;; dt =
;; rebin(reform(dt,1,1,dt_len),dem_map_sz[1],dem_map_sz[2],dt_len)

;multiply each intensity by dt so you only have to make one int map
;that is now scaled to be used with a DEM and not EM
for line_num = 0,n_elements(lines_small)-1 do begin
   lines_small[line_num].int *= dt

endfor



pad = fltarr(dem_map_sz[1],dem_map_sz[2]) ;to deal with periodicity of shift

TIC
;loop through each line and each dem basis element  
;; for i=0,basis_elements-1 do begin
   
;;    new_dem_map = meit_ll_map->scale_dem_map(dem_basis[*,i],/normalize)
   
;;    for j = 0,n_elements(lines_small)-1 do begin

;;       ;; int_map = rebin(reform(ll.lines[j].int,1,1,dt_len),dem_map_sz[1],dem_map_sz[2],dt_len)
;;       ;; disperse_pix = round(1/.028*(ll.lines(j).wvl-303.7860)) ; caluculate dispersion in pix
;;       int_map = rebin(reform(lines_small[j].int,1,1,dt_len),dem_map_sz[1],dem_map_sz[2],dt_len)
;;       disperse_pix = round(1/.028*(lines_small[j].wvl-303.7860)) ; caluculate dispersion in pix
 
;;       img = total(int_map*10.^(new_dem_map)*dt,3)
      
;;       ;scale image with MOSESI throughput curves and disperse
;;       zero = img*zero_throughput(j)
;;       plus = shift([img*plus_throughput(j),pad],[disperse_pix,0])
;;       minus = shift([img*minus_throughput(j),pad],[disperse_pix,0])

;;       ;assign to cube.
;;       moses_synth_cube[0,0,0,i] =  moses_synth_cube[*,*,0,i] + zero
;;       moses_synth_cube[0,0,1,i] =  moses_synth_cube[*,*,1,i] + plus[0:dem_map_sz[1]-1,*]
;;       moses_synth_cube[0,0,2,i] =  moses_synth_cube[*,*,2,i] + minus[0:dem_map_sz[1]-1,*]

;;       print,i,j
      
;;       TOC
    
;;    end
;; end




;when we use every point in temperature space, this calculation can be
;done a bit differently and is carried out here.


new_dem_map = meit_ll_map->get_dem_map()
basis_elements = dem_map_sz[3]

;find mean DEM to start
eit_dem = fltarr(basis_elements)
for i = 0,basis_elements-1 do eit_dem[i] = median(new_dem_map[0:2047,0:1023,i])

;normalize DEM by the average

new_dem_map -= rebin(reform(eit_dem,1,1,basis_elements),dem_map_sz[1],dem_map_sz[2],basis_elements) ;Subtract not Divide when working with Log DEM!!!

restore, 'mosesI_eit_alignment.sav'
for i=0,basis_elements-1 do begin
   
   for j = 0,n_elements(lines_small)-1 do begin

      disperse_pix = round(1/.028*(lines_small[j].wvl-303.7860)) ; caluculate dispersion in pix

      if lines_small[j].int[i] ne 0. then begin  ;avoid multiplying a bunch of zeros.
         img = lines_small[j].int[i]*10.^(new_dem_map[*,*,i])
         
                                ;scale image with MOSESI throughput curves and disperse
         zero = img*zero_throughput(j)
         plus = shift([img[0:-1-m[0]-1,*]*plus_throughput(j),pad,img[-1-m[0]:-1,*]*plus_throughput(j)],[-disperse_pix,0])  ;improved padding accounts for the periodicty of shift at the east limb
         minus = shift([img[0:-1-m[0]-1,*]*minus_throughput(j),pad,img[-1-m[0]:-1,*]*minus_throughput(j)],[disperse_pix,0])

                                ;assign to cube.
         moses_synth_cube[0,0,0,i] =  moses_synth_cube[*,*,0,i] + zero
         moses_synth_cube[0,0,1,i] =  moses_synth_cube[*,*,1,i] + plus[0:dem_map_sz[1]-1,*]
         moses_synth_cube[0,0,2,i] =  moses_synth_cube[*,*,2,i] + minus[0:dem_map_sz[1]-1,*]
      endif
      print,i,j
      
      TOC
    
   end
end


moses_synth_cube = moses_synth_cube[0:2047,0:1023,*,*]
restore, 'psfs2.sav'

  for i = 0,n_elements(moses_synth_cube[0,0,0,*])-1 do begin
     moses_synth_cube[0,0,0,i] = convol(moses_synth_cube[*,*,0,i],psfz)
     moses_synth_cube[0,0,1,i] = convol(moses_synth_cube[*,*,1,i],psfp)
     moses_synth_cube[0,0,2,i] = convol(moses_synth_cube[*,*,2,i],psfm)
  endfor

save, moses_synth_cube,eit_dem,filename="meit_synth_cube_full.sav"

end
