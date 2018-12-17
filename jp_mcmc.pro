;This will contain all the code neeeded to take the synthetic meit
;cube genereated by meit_synthetic and find the best fit to the mean
;MOSES difference image cross-correlation funtion

function meit_likelihood,param,data,test,tcor = tcor,dif_img=dif_img
;apply scaling of param to data
  param_sz = size(param)
  data_sz = size(data)

  data_temp = data
;to conserve memory, I am just going to do this in a loop so I
;don't have to form another massive array.
  for i = 0,param_sz[1]-1 do begin
     data_temp[*,*,*,i] *=10.^param[i]
  endfor


  data_temp = total(data_temp,4)           ;collapse the param dimension

;divide out the median of each image prior to subtraction

  for i = 0,data_sz[3]-1 do begin
     if  median(data_temp[*,*,i]) ne 0. then begin
        data_temp[*,*,i] /= median(data_temp[*,*,i])
      
     endif
     
  endfor
  
    

;compute cc
  img_sz = size(data[*,*,1])
  ;form difference images
  pz = data_temp[*,*,1]-data_temp[*,*,0]
  mz = data_temp[*,*,2]-data_temp[*,*,0]

  ;mean subtract
  pz -= make_array(img_sz(1),value=1)#mean(pz,dimension=1)
  mz -= make_array(img_sz(1),value=1)#mean(mz,dimension=1)

  if keyword_set(dif_img) then dif_img = pz
  ;pad with zeros
  pz_pad = [pz,fltarr(img_sz(1),img_sz(2))]
  mz_pad = [mz,fltarr(img_sz(1),img_sz(2))]


  f_pz = FFT(pz_pad,dimension=1)
  f_mz = FFT(mz_pad,dimension=1)
  test_cor = img_sz(1)*2*FFT(f_pz*conj(f_mz),dimension=1,/inverse)
  normalization = ((fltarr(n_elements(pz_pad(*,0)))+1)#sqrt(total(pz^2,1)*total(mz^2,1)))
  normalization[where(normalization eq 0)]=1e-6        ;avoid divide by zero error
  test_cor /= normalization ;normalize cc
  test_cor = reverse(real_part(test_cor))
  test_cor = shift(test_cor,[-img_sz(1)+1,0])
  tcor = mean(test_cor,dimension=2)
 

                                ;going to implement a real chi squared test for this
  normalization = total((tcor-mean(tcor))^2)
  error = total((test-tcor)^2)
  error = error/normalization
  
  return, error

end


;----------------------------------------------------------------------------
;Actual minimization routine, hoping to make this a little more
;generic.

; likelihood will be a string pointing to the likelihood
; funtion. Maybe inpliment this later

; param is the array of free parameters

; limits is an array the is n_elements(param)x 2 specifying some
; boundaries in parameter space.

; max_chain_length is when we cut off the Markov Chain if it
; hasn't converged (convergence criterion TBD so specfiy this
; so it stops)

; test is the function to compare the fit to

; if /random_start is set, initial param values will be ignored and
; replaced with random params

function jp_mcmc,param,data,limits,max_chain_length,test,width,random_start = random_start,verbose = verbose
TIC
 
  
  ;generate number between low and high limit
  if keyword_set(random_start) then param = randomu(seed,n_elements(param))*(limits[*,1]-limits[*,0])+limits[*,0]

  param_history = param         ;initialize keeping track of the parameters used
  
  error = meit_likelihood(param,data,test) ;initialize the first link in the chain
  print,error
;need to do a "burn in" I guess.
  ;; i = 0
  ;; burn_length = 100
  ;; Print,'Burning in!!'
  ;; while i lt burn_length do begin
  ;;    width = 1. 
  ;;    kick = randomu(seed,n_elements(param))*width-width/2 ;random kick between -1 and 1
  ;;    new_param = param+kick
  ;;    new_param = smooth(new_param,3)
  ;;    ;if kick puts you out of bounds, set it to the line
  ;;    where_less = where(new_param lt limits[*,0])
  ;;    where_more = where(new_param gt limits[*,1])
  ;;    if total(where_less) ne -1 then new_param[where_less]=limits[where_less] 
  ;;    if total(where_more) ne -1 then new_param[where_more]=limits[where_more]
  ;;    new_error = meit_likelihood(new_param,data,test)
  ;;    if new_error lt error[-1] then begin
  ;;       param_history = [[param_history],[new_param]]
  ;;       error = [error,new_error]
  ;;       param = new_param
  ;;       i+=1
  ;;    endif else begin
  ;;       error_ratio = (new_error-error[-1])/new_error
  ;;       does_it_stick = randomu(seed,1)
  ;;       if does_it_stick gt error_ratio then begin
  ;;          error = [error,new_error]
  ;;          param_history = [[param_history],[new_param]]
  ;;          param = new_param
  ;;          i += 1
           
  ;;       endif else print,"Parameters Rejected"
  ;;       ;If you get here that means new error was rejected
  ;;    endelse
        
        
  ;;    if keyword_set(verbose) then begin
  ;;       print,new_error
  ;;       print,burn_length-i
  ;;       TOC
  ;;    endif
  ;; endwhile

  
  ;; param_history = param_history[*,where(error eq min(error))]
  ;; param = param_history
  ;; error = error[where(error eq min(error))]
  
  
  i=0
  
  while i lt max_chain_length do begin
    
     ;; kick = randomu(seed,n_elements(param))*width-width/2
      kick = randomn(seed,n_elements(param))*width
     
     new_param = param+kick
     ;if kick puts you out of bounds, set it to the line
     where_less = where(new_param lt limits[*,0])
     where_more = where(new_param gt limits[*,1])
     if total(where_less) ne -1 then new_param[where_less]=limits[where_less] 
     if total(where_more) ne -1 then new_param[where_more]=limits[where_more]
     new_error = meit_likelihood(new_param,data,test)
     
     print,new_error
     if new_error lt error[-1] then begin
        error = [error,new_error]
        param_history = [[param_history],[new_param]]
        param = new_param
        i += 1
        
        endif else begin
        error_ratio = (new_error-error[-1])/new_error
        does_it_stick = randomu(seed,1)
        
        if does_it_stick gt error_ratio then begin
           error = [error,new_error]
           param_history = [[param_history],[new_param]]
           param = new_param
           i += 1
           
        endif else print,"Parameters Rejected"

     ;If you get here that means new error was rejected   
     endelse

     ;i+=1 ; do this for now so it stops running.
     

     if keyword_set(verbose) then begin
        print,max_chain_length - i
        TOC
     endif
     
  endwhile

  fit = {param_history:param_history,error_history:error}
  return,fit
  

end







