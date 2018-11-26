function linelistmap::init
  ;-- allocate memory to pointer when initializing object
  self.dem_map=ptr_new(/allocate)
  self.linelist=ptr_new(/allocate)
  return,1
end

;----------------------------------------------------------------

pro linelistmap::set_dem_map,value
;-- if data value exists, then insert into pointer location
  if n_elements(value) ne 0 then *(self.dem_map)=value
  return
end

;----------------------------------------------------------------

pro linelistmap::set_linelist,value
;-- if data value exists, then insert into pointer location
  if n_elements(value) ne 0 then *(self.linelist)=value
  return
end


;----------------------------------------------------------------

function linelistmap::get_linelist,value
;-- if data value is stored in object pointer, then copy it out
  if n_elements(*(self.linelist)) ne 0 then value=*(self.linelist)
  return,value
end

;------------------------------------------------------------------

function linelistmap::get_dem_map,value
;-- if data value is stored in object pointer, then copy it out
  if n_elements(*(self.dem_map)) ne 0 then value=*(self.dem_map)
  return,value
end
;------------------------------------------------------------------

function linelistmap::get_ion_image,snote
;-- build dt array for scaling from EM to DEM

;since ch_synthetic uses EM and not DEM in Isothermal mode we need to
;prepare an array of dt to multiply intensities before
;integrating. I don't know why there is a natural log, but it is
;copied from chianti ch_synthetic

  line_list = *(self.linelist)
  temp = 10.^line_list.logt_isothermal
  dlnt = ALOG(10.^(line_list.logt_isothermal[1]-line_list.logt_isothermal[0]))
  dt = temp * dlnt
 
;-- rebin and reform to multiply against dem map
  dem_map_sz = size(*(self.dem_map))
  dt_len = n_elements(dt)
  dt = rebin(reform(dt,1,1,dt_len),dem_map_sz[1],dem_map_sz[2],dt_len)

;-- find desired ion
  i  = where(line_list.lines.snote eq snote)

  if total(i) eq -1 then begin

     print,'No such ion in linelist!!!  Try again.'
  
     return,-1
  endif else begin
     if n_elements(i) gt 1 then begin
        int = total(line_list.lines[i].int,2)
        int_map = rebin(reform(int,1,1,dt_len),dem_map_sz[1],dem_map_sz[2],dt_len)
     endif else int_map = rebin(reform(line_list.lines[i].int,1,1,dt_len),dem_map_sz[1],dem_map_sz[2],dt_len)
  endelse

  image = total(int_map*10.^*(self.dem_map)*dt,3)
  return,image
end

function linelistmap::scale_dem_map,new_dem,normalize=normalize,overwrite=overwrite
;-- Take a new dem input and multiply it against the existing dem_map,
;   if keyword_set(normalize), divide the dem_map by its own median
;   DEM.  If overwrite set, set value of dem_map to new dem_map
  new_dem_map = *(self.dem_map)
  dem_map_sz = size(*(self.dem_map))
 
  if n_elements(new_dem) ne dem_map_sz[3] then begin
     print,'DEM does not have correct dimensions.'
     return,-1
  end

     
  
  if keyword_set(normalize) then begin
    
     median_dem = median(median(*(self.dem_map),dimension=1),dimension=1)
     new_dem_map -= rebin(reform(median_dem,1,1,dem_map_sz[3]),dem_map_sz[1],dem_map_sz[2],dem_map_sz[3])
  end

  new_dem_map -= rebin(reform(new_dem,1,1,dem_map_sz[3]),dem_map_sz[1],dem_map_sz[2],dem_map_sz[3])

  if keyword_set(overwrite) then self->set_dem_map,new_dem_map

  return,new_dem_map

end

  
  
pro linelistmap::cleanup
;-- free memory allocated to pointer when destroying object
  ptr_free,self.dem_map
  ptr_free,self.linelist
  return
end

;------------------------------------------------------------------

pro linelistmap__define

  void = {linelistmap, dem_map:ptr_new(), linelist:ptr_new()}
  return

end

  
