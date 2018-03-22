; Program to plot MOSES super exposures with a uniform scaling.  Used to make images for ParkerKankelborg.tex

restore,'moses_super.sav'



; Comment out old coyote graphics ... a thing of the past
  
;; cgwindow,'cgimage',super_plus,maxvalue=5,minvalue=0,ctindex=3,/keep_aspect_ratio,waspect=.5,wtitle='plus' 

;; cgwindow,'cgimage',super_zero,maxvalue=5,minvalue=0,ctindex=3,/keep_aspect_ratio,waspect=.5,wtitle='zero'

;;  eit_pz = super_plus-super_zero
;;  pz = (eit_pz-mean(eit_pz))
;;  pz /= sqrt(total(total(pz^2,2),1))

;; cgwindow,'cgimage',pz,maxvalue=.0005,minvalue=-.0005,ctindex=3,/keep_aspect_ratio,waspect=.5
;; cgwindow,'cgpolygon', [0.35, 0.35, 0.65, 0.65, .35], [0.7, 0.97, 0.97, 0.7, 0.7], /NORMAL, COLOR='green',thick='10',/addcmd
;; cgwindow,'cgpolygon', [0.65, 0.65, 0.95, 0.95, .65], [0.25, 0.5, 0.5, 0.25, 0.25], /NORMAL, COLOR='green',thick='10',/addcmd
;; cgwindow,'cgpolygon', [0.02, 0.02, 0.07, 0.07, .02], [0.02, 0.2, 0.2, 0.02, 0.02], /NORMAL, COLOR='green',thick='10',/addcmd,wtitle='Plus - Zero'

  plus = image(super_plus,max_value = 5, min_value = 0, rgb_table = 3,buffer=1,margin=0)
  zero = image(super_zero,max_value = 5, min_value = 0, rgb_table = 3,buffer=1,margin=0)

  eit_pz = super_plus-super_zero
  pz = (eit_pz-mean(eit_pz))
  pz /= sqrt(total(total(pz^2,2),1))
  dif_plot = image(pz,max_value =.0005, min_value = -0.0005, rgb_table = 3,buffer=1,margin=0)

 

  p = polygon([0.35, 0.35, 0.65, 0.65, .35], [0.65, 0.8, 0.8, 0.65, 0.65], /normal, target=dif_plot,fill_transparency=100,color='green',thick=3)
  p = polygon([0.65, 0.65, 0.95, 0.95, .65], [0.35, 0.5, 0.5, 0.35, 0.35], /normal, target=dif_plot,fill_transparency=100,color='green',thick=3)
  p = polygon([0.02, 0.02, 0.07, 0.07, .02], [0.2, 0.3, 0.3, 0.2, 0.2], /normal, target=dif_plot,fill_transparency=100,color='green',thick=3)
  

  plus.save, "../MCOR_paper/super_plus.eps", border=0
  zero.save, "../MCOR_paper/super_zero.eps", border=0
  dif_plot.save, "../MCOR_paper/super_pz.eps", border=0

  
end
