;This is my attempt to create red_noise of a given
;dimension...specified by dimension.  Alpha is the lag 1 auto
;correlation. White gives the white noise generated to make the red noise.

function red_noise,dimension,alpha,white=white


  
  white= randomu(seed,dimension,/normal)


  red =fltarr(n_elements(white))
  for n=1,n_elements(white)-1 do begin

     red(n)= alpha*red(n-1)+white(n)

  end

  return,red
  




end
