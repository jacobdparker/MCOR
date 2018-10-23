;build a basis of synthetic spectra for MCOR forward model

;build a list of .dem files

dems = file_search('./custom_dem/T*.dem')

;call chianti with specfici input parameters

wmin = 280
wmax = 340
density = 1e9
ioneq_name = '/usr/local/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'

for i = 0,n_elements(dems)-1 do begin
   

   ch_synthetic, wmin, wmax, density = density,ioneq_name = ioneq_name, dem_name = dems[i],$
                 output = line_list
   save_file = './custom_spectra/Triangle'+strtrim(strcompress(string(i)),1)+'.fits'
   ch_write_fits,line_list,save_file

endfor





end
