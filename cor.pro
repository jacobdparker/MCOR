pro cor,X,Y,Z

;Cross Correlation

fx = FFT(x)
cfy = conj(FFT(y))
cfy = fx*cfy
Z = float(fft(cfy,/inverse))
sz = size(x)
npix = N_elements(x)  ;For Some Reason Multiplying by the number of elements in the padded 	image normalizes all the FFT's and returns actual correlation

Z=npix*shift(Z,[sz(1)/2,sz(2)/2])   ;center the correlation so that it is easier to look at

end


