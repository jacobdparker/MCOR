;make a grid in log T from 4 to 7

logT = findgen(25)/(8)+4
logT_extent = logT[-1] - logT[0]


basis_elements = 7


;make a triangle function
width = logT_extent/(basis_elements-1)
r = width

for i = 0,basis_elements-1 do begin

   shift = logT[0]+r*i
  
   y = ( (1-abs(float(logT-shift)/r)) > 0 )
   file_name = './custom_dem/Triangle'+strtrim(strcompress(string(i)),1)+'.dem'
   write_csv,file_name,logT,y

   ;append a "referenence" to each file as is expected by chianti
   spawn0 = 'echo "-1" >> '+file_name
   spawn1 = 'echo "This dem was created custom by Jake Parker for a MOSES forward model" >> '+file_name
   spawn2 = 'echo "-1" >> '+file_name
   spawn, spawn0
   spawn, spawn1
   spawn, spawn2

   if i eq 0 then plot,logT,y else oplot,logT,y

endfor








end
