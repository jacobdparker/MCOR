min_int = .001

  f = file_search('ch_sp/*_sp*')

  images = moses1_eit_ch(min_int,file_search('ch_sp/*ar.*'),ions=ions,ar_wvl=wvl,ar_int=int)
  ar = {mag:0.0,images:images,ions:ions,wvl:wvl,int:int}
  
images = moses1_eit_ch(min_int,file_search('ch_sp/*qs.*'),ions=ions,ar_wvl=wvl,ar_int=int)
  qs =  {mag:0.0,images:images,ions:ions,wvl:wvl,int:int}

  images = moses1_eit_ch(min_int,file_search('ch_sp/*qs_eis.*'),ions=ions,ar_wvl=wvl,ar_int=int)
  qs_eis = {mag:0.0,images:images,ions:ions,wvl:wvl,int:int}

  images = moses1_eit_ch(min_int,file_search('ch_sp/*chole.*'),ions=ions,ar_wvl=wvl,ar_int=int)
  chole =  {mag:0.0,images:images,ions:ions,wvl:wvl,int:int}

  images = moses1_eit_ch(min_int,file_search('ch_sp/*prom.*'),ions=ions,ar_wvl=wvl,ar_int=int)
  prom =  {mag:0.0,images:images,ions:ions,wvl:wvl,int:int}

  images = moses1_eit_ch(min_int,file_search('ch_sp/*flare.*'),ions=ions,ar_wvl=wvl,ar_int=int)
  flare =  {mag:0.0,images:images,ions:ions,wvl:wvl,int:int}
  
images = moses1_eit_ch(1000000,file_search('ch_sp/*qs.*'),ions=ions,ar_wvl=wvl,ar_int=int)
  he_ii =  {mag:0.0,images:images,ions:ions,wvl:wvl,int:int}
  

eit_ch = {he_ii:he_ii,qs:qs,qs_eis:qs_eis,ar:ar,flare:flare,prom:prom,chole:chole}




save,eit_ch,filename='eit_ch_all.sav'

end
