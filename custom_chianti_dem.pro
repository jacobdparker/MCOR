; A few notes.
; ref = string describing the dem, could probably just input a blank
; string

; dem = chianti expects these to actually be log dem in units cm^-5 K^-1

pro custom_chianti_dem,dem_logt,dem,ref,file_name = file_name


  
  if keyword_set(file_name) eq 0 then begin
     file_name = ''
     read, file_name, PROMPT='Enter path to new .dem: '

     ;add possible error checking if file doesn't end in .dem
  end

   write_csv,file_name,dem_logT,dem

   ;append a "referenence" to each file as is expected by chianti
   spawn0 = 'echo "-1" >> '+file_name
   spawn1 = 'echo '+ref+'>>'+file_name
   spawn2 = 'echo "-1" >> '+file_name
   spawn, spawn0
   spawn, spawn1
   spawn, spawn2

   

end
