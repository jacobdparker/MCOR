

;; ;Main driver fo the minimization scheme.  Things to keep in mind,
;; ;right now parameters are in log, to better match the range of values
;; ;we expect in DEMs  Should not exceed 10 I don't think.


;;   restore, 'meit_synth_cube.sav'
;;   restore, 'mcor.sav'
  
;;   param = fltarr(n_elements(moses_synth_cube[0,0,0,*]))
;;   limits = [[param],[param]]
;;   limits[*,0]=-10                ;ten orders of magnitude difference between dem basis
;;   max_chain_length = 100

;;   data = moses_synth_cube
;;   test = mcor
;;   restore,'mcmc_goodfit.sav'
  
;;   fit = jp_mcmc(param,data,limits,max_chain_length,test,/verbose)

;;   tcor = mcor                   ;start with something
;;   dif_img = 1
;;   param = fit.param_history[*,where(fit.error_history eq min(fit.error_history))]
;;   error = meit_likelihood(param,data,test,tcor=tcor,dif_img = dif_img)

  ;recreate DEM basis, for now triangle functions in logt/logdem space
  restore, 'moses_eit_linelist_map.sav'
  ll = meit_ll_map->get_linelist()
  logt = ll.logt_isothermal
  logt_sz = n_elements(logt)

  basis_elements = n_elements(param)
  height = 10
  width = (logt_sz-1)/(basis_elements-1)
  
  x = findgen(logt_sz)
 


  dem_basis = fltarr(logt_sz,basis_elements)
  for i = 0,basis_elements-1 do begin
     dem_basis[0,i] = height-abs((x-width*i)*height/width)
     dem_basis[where(dem_basis[*,i] le 0),i] = 0
  end
 
 
 
  
  STOP
  dem_basis += rebin(reform(param,1,n_elements(param)),n_elements(logt),n_elements(param))
  best_dem = total(dem_basis,2)

end





