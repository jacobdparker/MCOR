function meit_error,param,tcor = tcor,dif_img=dif_img
    common amoeba, data, test
    
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

    alpha = 1
    smooth_kernel = 5
    error = total(abs(tcor-test)) + alpha * total(abs(smooth(param,smooth_kernel)-param))
    ;/total(abs(tcor))


    return, error
end

pro meit_synthetic_amoeba
    common amoeba, data, test

    TIC

    restore, 'meit_synth_cube_full.sav'
    restore, 'mcor.sav'
 
    data = moses_synth_cube
    test = mcor
    param = eit_dem
 
    


    ;some manual guessing to get us started
    param[0:6] += 3.55
    param[20] += .75
    param[15] -= 2.5

    limits = [[param],[param]]
    limits[*,0]=15               
    limits[*,1]=30

    ;random start
    param = randomu(seed,n_elements(param))*(limits[*,1]-limits[*,0])+limits[*,0]

    

    ;amoeba input parameters
    ftol = 1.   
    fname = 'meit_error'
    scl = 3
    nmax = 5000
    
    ;find minima from eit_dem
    best_fit = amoeba(ftol,function_name = fname , function_value = error, scale = scl, p0 = param, nmax = nmax)
    best_fit_error = error[0]
    print,error[0]
    

    ;relaunch amoeba from end point
    flag = 0
    while flag eq 0 do begin
        result = amoeba(ftol,function_name = fname , function_value = error, scale = scl, p0 = best_fit, nmax = nmax)
        if error[0] lt best_fit_error then begin
            best_fit_error = error[0]
            best_fit = result
            print, best_fit_error
        endif else begin
            flag = 1 
            print, "No better fit found"
            print, "Best Fit Error = ",best_fit_error
            plot, best_fit
        endelse

    end
         
        
        
        


    ;n_amoeba = 10
    ;;launch a bunch of random starting points
    ;for i = 0,n_amoeba do begin
    ;    param = randomu(seed,n_elements(param))*(limits[*,1]-limits[*,0])+limits[*,0]
    ;    result = amoeba(ftol,function_name = fname , function_value = error, scale = scl, p0 = param, nmax = nmax)
    ;    best_fit_error = [best_fit_error, error[0]]
    ;    best_fit = [[best_fit],[result]]
    ;    print, error[0]
    ;end

    TOC
    STOP
end

    


