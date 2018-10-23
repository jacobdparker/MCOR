; build an eit_line map for the MOSES I launch

f171 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.190014'
f195 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.183429'
f284 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.182651'
f304 = '/home/jake/Documents/MOSES/MCOR/moses_eit/efz20060208.184020'


restore, 'eit_cube_prepped.sav'

eit_dem_tool,eit_cube,dem_map,temp,/mk_line,ion = 'si11', line_map = line_map



;eit_dem_tool,[f171,f195,f284,f304],dem_map,temp,/mk_line, $
            ; ion = 'si11', line_map = line_map

end
