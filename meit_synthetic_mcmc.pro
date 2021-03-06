

;Main driver for the minimization scheme.

;Instructions
; .compile jp_mcmc
; .r meit_synthetic_mcmc
TIC

  restore, 'meit_synth_cube_full.sav' ;created by meit_synthetic
  restore, 'mcor.sav'
 

  param = eit_dem
 
  ;set the bounds for each bin in LOG DEM space
  limits = [[param],[param]]
  limits[*,0]=15               
  limits[*,1]=30

  ;; fixd = 0
  ;; limits[0:fixd,0] = param[0:fixd]
  ;; limits[0:fixd,1] = param[0:fixd]


  ;some manual guessing to get us started
  param[0:6] += 3.55

  param[20] += .75

  param[15] -= 2.5


 
 
  max_chain_length = 10000
  width = .5

  
  data = moses_synth_cube
  test = mcor
  STOP
                            




  cd, current=dir

  n_proc = 4

  if n_proc gt 1 then begin
     
  
     br2 = objarr(n_proc)
     for i=0, n_proc-1 do begin
                                ;needed for any code like this
        br2[i] = IDL_IDLbridge()
        br2[i]->SetVar,'dir',dir
        br2[i]->Execute,'cd ,dir'
        br2[i]->SetVar,'!path',!path

                                ;MCMC specficic initialize inputs
        br2[i]->SetVar,'param',param
        br2[i]->SetVar,'data',data
        br2[i]->SetVar,'limits',limits
        br2[i]->SetVar,'max_chain_length',max_chain_length
        br2[i]->SetVar,'test',test
        br2[i]->SetVar,'width',width
        
                                ;Run Chain
        br2[i]->Execute,'fit = jp_mcmc(param,data,limits,max_chain_length,test,width,/verbose,/random_start)',/nowait
     endfor

                                ;Monitor Chain Status
     status = replicate(1, n_proc)
     while max(status) gt 0 do begin
        wait, 1
        for i=0, n_proc-1 do begin
           status[i] = br2[i]->Status()
        endfor
        
     endwhile

     br2[0]->Execute,'error_history = fit.error_history'
     br2[0]->Execute,'param_history = fit.param_history'
     fit = {error_history:br2[0]->GetVar('error_history'),param_history:br2[0]->GetVar('param_history')}
     fits = fit
     for i =1,n_proc-1 do begin
        br2[i]->Execute,'error_history = fit.error_history'
        br2[i]->Execute,'param_history = fit.param_history'
        fit = {error_history:br2[i]->GetVar('error_history'),param_history:br2[i]->GetVar('param_history')}
        fits = [fits,fit]
     endfor
     
     for j=n_proc-1,0,-1 do obj_destroy, br2[j]

     min_errors = fltarr(n_proc)
     for i = 0,n_proc-1 do min_errors[i]=min(fits[i].error_history)
     lowest_error = min(min_errors)
     best_chain = where(min_errors eq lowest_error)

     if n_elements(best_chain) gt 1 then fit = fits[0] else fit = fits[best_chain]

    
     
  endif else begin
     
     fit = jp_mcmc(param,data,limits,max_chain_length,test,width,/verbose,/random_start)
     lowest_error = min(fit.error_history)
  endelse 
  
  

  
  tcor = mcor                   ;start with something for tcor so we can get it back out
  dif_img = 1
  param = fit.param_history[*,where(fit.error_history eq lowest_error)]
  error = meit_likelihood(param,data,test,tcor=tcor,dif_img = dif_img)

  restore,'moses_eit_linelist_map.sav'
  ll = meit_ll_map->get_linelist()

  DEM_plot = plot(ll.logt_isothermal,param,DIMENSIONS=[1030,800])
  DEM_plot.title = 'Mean EIT DEM'
  DEM_plot.ytitle = '$log_{10}$ DEM (arbitrary units)'
  DEM_plot.xtitle  = '$log_{10}$ Temperature ($log_{10}$ K)'
  DEM_plot.thick = 4
  DEM_plot.font_size = 20





  print,'Lowest Error  =', error

  restore,'mosesI_pointing.sav' 
  
  dif_img= (dif_img-mean(dif_img))
  dif_img /= sqrt(total(total(dif_img^2,2),1))
  dif_plot = image(dif_img,x,y,max_value =.0005, min_value = -0.0005, rgb_table = 3,buffer=1,axis_style=1,title = 'Modeled MOSES Difference Image (m=+1 minus m=0)', xtitle ='X (arcsec)',ytitle = 'Y (arcsec)')

  lag = indgen(n_elements(mcor))
  lag -= max(lag)/2
  
  
  cc_plot = plot(lag,mcor,name = 'MOSES',color = 'b')
  cc2_plot = plot(lag,tcor,color='r',/overplot,name = 'Modeled')
  leg = LEGEND(TARGET=[cc_plot,cc2_plot],position = [-500,.8], $
  /DATA, /AUTO_TEXT_COLOR)

  cc_plot.title = 'Difference Image Cross-Correlation'
  cc_plot.xtitle = 'Lag (pix)'
  cc_plot.ytitle = 'Correlation'
  cc_plot.xrange = [-1950,1950]

TOC
  
  ;save plots
  ;; dif_plot.save, "../AGU_2018/figures/modeled_pz.eps"
  ;; cc_plot.save,  "../AGU_2018/figures/cc_comp.eps"
  ;; dem_plot.save, "../AGU_2018/figures/eit_dem.eps"
end





