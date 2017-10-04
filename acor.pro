function acor,X,Y,dimension=dimension



;Auto Correlation
fx = FFT(x,dimension=dimension)
cfx =conj(fx)
cfx= fx*cfx
cfx=FFT(cfx,dimension=dimension,/INVERSE)
Y=float(cfx)
sz = size(x)
npix = N_elements(x)  ;For Some Reason Multiplying by the number of elements in the padded 	image normalizes all the FFT's and returns actual correlation

Y=npix*shift(Y,[sz(1)/2,sz(2)/2])   ;center the correlation so that it is easier to look at

Y = Y / total(Y^2,dimension)


return, Y
end

