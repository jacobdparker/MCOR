

PRO throughput,stack,y,t,h

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  NOTE!!! Possibly a flaw in how I am dividing out the "overlap"
;  figure that shit out
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Subtract mean of Image to make auto correlation viewable

stack=stack-mean(stack)

;Zero Pad and Autocorrelate Image

pad,stack,x
acor,x,y


;make a matrix with entries of the 'overlap'
sz=size(stack)
t=[findgen(sz(2)),rotate(findgen(sz(2)+1),2)]
t=float(transpose(t))
tt=[findgen(sz(1)),rotate(findgen(sz(1)+1),2)]
t=t##tt

;
;Find Grid Height
;


h=fltarr(64)  ;build array to contain grid heights

sz=size(y)
center=[sz(1)/2,sz(2)/2]  ;start in the center

;build array containing various maximum positions
mx=fltarr(2,1,64)

blah=0
for m=0,7 do begin
for q=1,8 do begin
mx[0,0,blah]=center+q*[27,-4]+m*[3,27]  ;fill in location of various maxima
blah=blah+1
end 
end

print,'max coordinates'
print,mx

for l=0,63 do begin

;Build a mean MAX value by grabbing 9 pixels around apparent max
m=fltarr(3,3)
for i=-1,1 do begin     
for j=-1,1 do begin

m[i+1,j+1]=y[mx[0,0,l]+i,mx[1,0,l]+j]/t[mx[0,0,l]+i,mx[1,0,l]+j]

end
end
;print,'max matrix'
;print,m

mx1=mean(m)
;print,'Max'
;print,mx1


;Begin gathering surrounding minima
minima=fltarr(4)
pos=fltarr(2,1,4)  ;array of surrounding minima positions
pos[0,0,0]=mx[0:1,0,l]+[7,6]
pos[0,0,1]=mx[0:1,0,l]+[-7,-6]
pos[0,0,2]=mx[0:1,0,l]+[6,-7]
pos[0,0,3]=mx[0:1,0,l]+[-6,7]

for n=0,3 do begin  
m=fltarr(3,3)
for i=-1,1 do begin     
for j=-1,1 do begin

m[i+1,j+1]=y[pos[0,0,n]+i,pos[1,0,n]+j]/t[pos[0,0,n]+i,pos[1,0,n]]

end
end
minima[n]=mean(m)
end
mn1=mean(minima)
;print,'min'
;print,mn1

h[l]=mx1-mn1

end


h=h/.82
h=sqrt(h)

print,'GRID HEIGHT'
print,h

end
