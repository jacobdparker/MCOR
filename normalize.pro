function normalize,x

sz = size(x)
y = fltarr(sz(1),sz(2),sz(3))

for n = 0,sz(3)-1 do begin
y(*,*,n)= x(*,*,n)/median(x(*,*,n))
end 

return,y

end
