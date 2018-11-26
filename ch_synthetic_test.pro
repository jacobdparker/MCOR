; test differences between isothermal and non-isothermal treatments in
; ch_synthetic


wmin = 280
wmax = 340
pressure = 1e15
ioneq_name = '/usr/local/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'

read_ioneq,ioneq_name,ioneq_logt,ioneq,ioneq_ref ;read in ioneq_logt to avoid interpolation

dem_name = !xuvtop+'/dem/quiet_sun.dem'
read_dem,dem_name,dem_logt,dem,ref

temp = 10.^dem_logt
dlnt=ALOG(10.^(dem_logt[1]-dem_logt[0]))
dt=temp * dlnt





ch_synthetic, wmin, wmax, pressure=pressure,dem_name = dem_name, ioneq_name = ioneq_name,$
              output = line_list_dem,sngl_ion='He_2'


ch_synthetic, wmin, wmax, pressure=pressure,logt_isothermal = dem_logt, ioneq_name = ioneq_name,$
              output = line_list_isothermal,sngl_ion='He_2',/no_sum_int




end
