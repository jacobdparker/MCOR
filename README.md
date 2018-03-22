# MCOR
Software for identifying extra spectral content in the MOSES I Data

Since I did a shitty job of documenting my work I get to do it retroactively.  My goal is to identify the final versions of various pieces of code, and the order they were run in.  Hopefully I can also delete old unused versions of the code so this package is useful to the MOSES team in the future.

#Forming Time Averaged Difference Images
The IDL routine 'moses_super_exporsure.pro' requires mosesAlignFinal.sav (a coaligned and reduced version of the MOSESI dataset).  Script output is moses_super.sav which contains super_zero, super_minus, and super_plus (a single image in each order with no saturated pixels.  More documentation within.

Plots of these super exposures are generated through 'moses_super_plot.pro'.  This saves the plus, zero, and plus-zero images with a uniform image scaling and draws some green boxes around interesting areas.

#Cross Correlation
... this part is extra bad.  mcor.sav is the cross-correlation function that I used for all forward modeling.  As far as I can tell sigtest_plot.pro will tell me everything I need to know if I can just parse through it.




