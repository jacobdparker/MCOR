restore,'moses_super.sav'

cgwindow,'cgimage',super_plus,maxvalue=5,minvalue=0,ctindex=3,/keep_aspect_ratio,waspect=.5,wtitle='plus'

cgwindow,'cgimage',super_zero,maxvalue=5,minvalue=0,ctindex=3,/keep_aspect_ratio,waspect=.5,wtitle='zero'

 eit_pz = super_plus-super_zero
pz = (eit_pz-mean(eit_pz))
     pz /= sqrt(total(total(pz^2,2),1))

cgwindow,'cgimage',pz,maxvalue=.0005,minvalue=-.0005,ctindex=3,/keep_aspect_ratio,waspect=.5
    cgwindow,'cgpolygon', [0.35, 0.35, 0.65, 0.65, .35], [0.7, 0.97, 0.97, 0.7, 0.7], /NORMAL, COLOR='green',thick='10',/addcmd
cgwindow,'cgpolygon', [0.65, 0.65, 0.95, 0.95, .65], [0.25, 0.5, 0.5, 0.25, 0.25], /NORMAL, COLOR='green',thick='10',/addcmd
cgwindow,'cgpolygon', [0.02, 0.02, 0.07, 0.07, .02], [0.02, 0.2, 0.2, 0.02, 0.02], /NORMAL, COLOR='green',thick='10',/addcmd,wtitle='Plus - Zero'

end
