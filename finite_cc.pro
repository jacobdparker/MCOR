;; compute the finite cross correlation of two vectors for given lag
;; in a way that I understand and control

function  finite_cc,x,y,lag

;pad arrays
sz = size(x)
xpad = [x,fltarr(sz(1))]
ypad = [y,fltarr(sz(1))]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Uncomment if you want a constant mean subtracted
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

xbar=mean(x)
ybar=mean(y)


;;;;;;;;
;begin cross correlation
;;;;;;;;;

lagsize = size(lag)



cc = fltarr(lagsize(1))
for n=0,lagsize(1)-1 do begin

yshift=shift(ypad,lag(n))

;;To eliminate the problem of negative indexing if statement changes
;;correltation for negative or postive lag

if lag(n) GE 0 then begin       ; section used for postive lag
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;this formulation has a mean that changes for a given lag;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; ;calculate mean of sections being correlated

;; xbar=mean(xpad(lag(n):sz(1)-1))
;; ybar=mean(yshift(lag(n):sz(1)-1))

;subtract that from the section being correlated
normx=xpad(lag(n):sz(1)-1)-xbar
normy=yshift(lag(n):sz(1)-1)-ybar


;; cc(n)=total(normx*normy)/sqrt(total(normx^2)*total(normy^2))

;;;;;;;;;;;;;;;;
;used for constant normalization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cc(n)=total(normx*normy)/sqrt(total((x-xbar)^2)*total((y-ybar)^2))


endif else begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;this formulation has a mean that changes for a given lag;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
;; ;calculate mean of sections being correlated

;; xbar=mean(xpad(0:sz(1)-1+lag(n)))
;; ybar=mean(yshift(0:sz(1)-1+lag(n)))



;subtract that from the section being correlated
normx=xpad(0:sz(1)-1+lag(n))-xbar
normy=yshift(0:sz(1)-1+lag(n))-ybar

;; cc(n)=total(normx*normy)/sqrt(total(normx^2)*total(normy^2)) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;used for constant normalization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
cc(n)=total(normx*normy)/sqrt(total((x-xbar)^2)*total((y-ybar)^2))

endelse




endfor

return,cc

end

