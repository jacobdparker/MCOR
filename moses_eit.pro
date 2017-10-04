;run moses_super_exposure

; Read in all of the eit image from the MOSES launch
eit_prep,'/disk/data/jparker/moses_eit/efz20060208.130013',index,i171
read_eit,'/disk/data/jparker/moses_eit/efz20060208.130013',i171_h

eit_prep,'/disk/data/jparker/moses_eit/efz20060208.130609',index,i284
read_eit,'/disk/data/jparker/moses_eit/efz20060208.130609',i284_h

eit_prep,'/disk/data/jparker/moses_eit/efz20060208.131348',index,i195
read_eit,'/disk/data/jparker/moses_eit/efz20060208.131348',i195_h

eit_prep,'/disk/data/jparker/moses_eit/efz20060208.131938',index,i304
read_eit,'/disk/data/jparker/moses_eit/efz20060208.131938',i304_h


i304 = congrid(i304,4512,4512)
i304 = i304(0:4095,0:4095)

;pad super_zero
pad_z = fltarr(4096,4096)
pad_z(0,0) = super_zero

;correlate both images
cor = convol_fft(i304,pad_z,/correlate)
cor = shift(cor,[-2047,-2047])

m = array_indices(pad_z,where(cor eq max(cor)))

i304= shift(i304,[-m(0),-m(1)])
i304=i304(*,0:1023)

i171 = congrid(i171,4512,4512)
i171 = i171(0:4095,0:4095)

i284 = congrid(i284,4512,4512)
i284 = i284(0:4095,0:4095)

i195 = congrid(i195,4512,4512)
i195 = i195(0:4095,0:4095)

i171= shift(i171,[-m(0),-m(1)])
i171=i171(*,0:1023)

i195= shift(i195,[-m(0),-m(1)])
i195=i195(*,0:1023)

i284= shift(i284,[-m(0),-m(1)])
i284=i284(*,0:1023)
end
