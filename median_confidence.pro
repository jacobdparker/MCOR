

;this program will calculate non tabulated (and tabulated actually) values for distribution
;free median confidence intervals based on equations in paper by
;J. L. Van der Parren



; n, the total sample size
; Confidence, defines the confidence interval
;
; Calling Sequence
; k=median_confidence(n,alpha,[Confidence=Confidence])
; k, index of data points that contain the confidence interval




function median_confidence,n,alpha,Confidence=Confidence,I=I


 
     
p = .5 ;this program I believe is written in a way that the p quantile could be given as an input.  This would require understanding what this means from the paper stated above.
k=0




if n le 170 then begin

   r =findgen(n-k+1)+k
   choose_vector =factorial(n)/factorial(n-r)/factorial(r)*p^(r)*(1-p)^(n-r)

   
   repeat begin

     k=k+1
     choose = total(choose_vector(k:n_elements(choose_vector)-1))
     
   endrep until choose le 1-alpha

   k=k-1
   choose = total(choose_vector(k:n_elements(choose_vector)-1))
   
endif else begin

   print,'!!! WARNING: Using Stirlings Approximation since n exceeds 170 !!!!'
   k=1 ;k=0 will throw an error with alog so we start with a higher guess
   n=double(n)
   
   r =double(findgen(n-k)+k)            ;here I have elimintated the n=k term in the sum to prevent alog(0) I think it is zero anyway
      
   ;take the ln for working with large numbers
   choose_vector = n*alog(n)-n + .5*alog(2*!pi*n)-(r*alog(r)-r + .5*alog(2*!pi*r))-((n-r)*alog(n-r)-(n-r) + .5*alog(2*!pi*(n-r))) + r*alog(p) + (n-r)*alog(1-p)
   choose_vector = exp(choose_vector)

   
   
   j=0ULL
   repeat begin
                                
      j=j+1
      choose = total(choose_vector(j:n_elements(choose_vector)-1))

      
   endrep until choose le 1-alpha


   j=j-1
   k=k+j
   choose = total(choose_vector(j:n_elements(choose_vector)-1))
endelse



I = 1-choose
Confidence = 1 - 2*I

return,k


end
