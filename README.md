# MCOR
Software for identifying extra spectral content in the MOSES I Data

Since I did a shitty job of documenting my work I get to do it retroactively.  My goal is to identify the final versions of various pieces of code, and the order they were run in.  Hopefully I can also delete old unused versions of the code so this package is useful to the MOSES team in the future.

#Forming Time Averaged Difference Images
The IDL routine 'moses_super_exporsure.pro' requires mosesAlignFinal.sav (a coaligned and reduced version of the MOSESI dataset).  Script output is moses_super.sav which contains super_zero, super_minus, and super_plus (a single image in each order with no saturated pixels.  More documentation within.

Plots of these super exposures are generated through 'moses_super_plot.pro'.  This saves the plus, zero, and plus-zero images with a uniform image scaling and draws some green boxes around interesting areas. Note:  I beilieve that image labeling is now all done in Latex, not by this script.

#Co-aligned EIT images
All done by running moses_eit.pro.  This takes the 4 full disk EUV images take by EIT close to launch and then despikes, rotates, and rebins to match MOSES FOV and resolution.  Final versions then saved in moses_eit.sav.

#Cross Correlation
... this part is really bad.  mcor.sav is the cross-correlation function that I used for all forward modeling.  As far as I can tell sigtest_plot.pro will tell me everything I need to know if I can just parse through it.

#Forward Model
Four prepped eit images described above are used with the eit_dem_tool to get a dem for every pixel.
An isothermal linelist is generated using the ch_synthetic routine from Chianti giving an intensity at every temperature bin for each line in the MOSES pass band.  
This is all controlled by meit_eit_linelist_map.pro.
A synthetic MOSES image for every temperature bin is created via meit_synthetic.pro .
Using the eit dem map, and an iso thermal linelist, a single line image for every temperature bin is shifted according to moses dispersion and then summed together.
This takes us from an image for every spectral line, at each temperature bin, to a single image per bin.  
Each bin is weighted by "dt" and the entire cube is divided by the median eit_dem, that way a synthetic MOSES image can be generated from any DEM by weighting the cube by a chose DEM and summing in temperature.
From there synthetic difference images can be made and cross correlated for comparison to MOSES data.


#Minimization
Right now my minimization method (mcmc) and merit function both exisit in jp_mcmc.pro
The driver program is meit_synthetic_mcmc.pro . 
The driver program sets parameter limits, determines the number of chains to run in parralel, and also makes some plots.






