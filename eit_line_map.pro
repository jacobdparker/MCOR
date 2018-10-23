;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;; NOTE: This file has been modified to point to the  appropriate  
;;;;;;; Chianti Files                                             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
; NAME:
;       EIT_LINE_MAP
;
; PURPOSE:
;       Create a line/bandpass map based upon input differential emission
;       measure (DEM) for EIT like instrument.
;
; CATEGORY:
;       Analysis
;
; CALLING SEQUENCE:
;       line_map =eit_line_map(wave,dem,intemp,x,coefs,ion=ion,wmin=wmin,$
;                      wmax=wmax,he_fac=he_fac,trange=trange,instr=instr,$
;                      irradiance=irradiance,delta_wave=delta_wave, $ 
;   			     xyidx=xyidx,nocon=nocon,abund=abund, $ 
;			     eqion=eqion,pressure=pressure)
;
; INPUTS:
;       wave                 - EIT bandpass or particular line wavelength
;       dem                  - input DEM (1D or 3D)
;       intemp              - temperature array (1D) corresponding to DEM
;       x                       - EIT specific line parameter
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       coefs                 - coefficients to "tweak" input DEM
;       ion                    - choice of output line if using single line
;       wmin                - minimum wavelength to consider for line maps
;       wmax                - maximum wavelength to consider for line maps
;       he_fac               - enhance factor for Helium II 304 over default
;                                  spectrum
;       trange                - array of temperature indices to consider if not
;                                   full range
;       instr                   - set to particular instrument, default = 'eit'
;                                   set to 'dum' or anything else for simple lines
;       irradiance          - set if only want line irradiance
;       delta_wave        - wavelength bin size for subsetting wavelength
;                                    range. Default is 1 bin for wmin->wmax. Only
;                                    useable with irradiance keyword
;       xyidx                 - Index values of array that should be used. All
;                                   others are set to zero
;       nocon                - Set if do not want continuum included in calculation
;       pressure             - Set if want spectra calculated with specific Pressure,
;                                   default = 1.0E+15
;       abund                - Set if want spectra calculated with specific Abundance,
;                                   default = Coronal (Feldman) ext
;       eqion                 - Set if want spectra calcualted with specific ionization
;                                   equilibrium, default = Mazzotta et al. ext
;
; OUTPUTS:
;       line_map            - computed line/bandpass map
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;
; COMMON BLOCKS: none.
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
; Requires Pre-computed spectra using CHIANTI Isothermal procedure
;
; PROCEDURE:
;
; SUBROUTINES: Include MK_EIT_SPEC for computing Isothermal spectra
;                              
; Helium enhancement = 9.5    Any change in the elemental abundances
;                                                employed from our Feldman default would
;                                                require re-examination of the Helium
;                                                enhancement factor
;
;
; MODIFICATION HISTORY:
;   Written by: J. Newmark  October 2001
;   Modified:	J. Newmark  May 2006
;   Modified:	J. Newmark  July 2010
;   Modified:	J. Cook	    Sept 2011
;   Modified:   J. Parker   October 2018: Edited default ioneq file
;-----------------------------------------------------------
pro mk_eit_spec,cont=cont,he_fac=he_fac,abund=abund,eqion=eqion,pressure=pressure,file=file,ion=ion
; Use CHIANTI Isothermal routine for calculating Unit EM spectra for each
;   temperature bin (4.0-6.5), optionally include continuum
;
; Atomic physics defaults
if n_elements(abund) eq 0 then $
    abund = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_ext.abund')
if n_elements(eqion) eq 0 then $
   eqion = concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq') 
   ;eqion = concat_dir(concat_dir(!xuvtop,'ioneq'),'mazzotta_etal_ext.ioneq')
if n_elements(pressure) eq 0 then pressure = 1.e15
if n_elements(he_fac) eq 0 then he_fac = 9.5
if keyword_set(cont) then begin
    cont = 1
    app = 'con'
endif else begin
    cont= 0
    app = 'nocon'
endelse
temp=indgen(26)*0.1+4.0
temp=10^temp

; For making dem_map, need isothermal with photons units
isothermal,165,550,1.,temp,waves,spec0,pressure=pressure,/noverbose,cont=cont,$
    abund=abund,ioneq=eqion
val304=where(waves eq 304)
spec0(val304,*)=spec0(val304,*)*he_fac
if n_elements(file) eq 0 then file = 'radiance_unit_dem_'+app+'.save'
save,file=file,temp,pressure,spec0,waves

if keyword_set(ion) then begin
; Here are three sample line calculations for He II 304 A, Si XI 303 A, Fe XIV 5303 A
; If you want another line, add new code appropriate for new line; keyword sngl in
;   isothermal call of format element_ionization stage
    isothermal,303.,304.,.02,temp,waves,spec0,pressure=pressure,/noverbose,cont=cont,$
        sngl='he_2',abund=abund,ioneq=eqion,/ergs
    save,file='he2_radiance_unit_dem_'+app+'.save',temp,pressure,spec0,waves

    isothermal,303.,304.,.02,temp,waves,spec0,pressure=pressure,/noverbose,cont=cont,$
        sngl='si_11',abund=abund,ioneq=eqion,/ergs
    save,file='si11_radiance_unit_dem_'+app+'.save',temp,pressure,spec0,waves

   isothermal,5304,5305,.02,temp,waves,spec0,pressure=pressure,/noverbose,cont=cont,$
       sngl='fe_14',abund=abund,ioneq=eqion,/ergs
   save,file='fe14_radiance_unit_dem_'+app+'.save',temp,pressure,spec0,waves

; To calculate line_map for broadband interval over 165-550 A range,
; want isothermal with ergs units
   isothermal,165,550,1.,temp,waves,spec0,pressure=pressure,/noverbose,cont=cont,$
    abund=abund,ioneq=eqion,/ergs
   val304=where(waves eq 304)
   spec0(val304,*)=spec0(val304,*)*he_fac
   save,file='broadband_radiance_unit_dem_'+app+'.save',pressure,spec0,waves
endif

end
;

function eit_line_map,wave,dem,intemp,x,coefs,ion=ion,wmin=wmin,wmax=wmax,$
          he_fac=he_fac,trange=trange,instr=instr,irradiance=irradiance,$
          delta_wave=delta_wave,xyidx=xyidx,nocon=nocon,$
          abund=abund,eqion=eqion,pressure=pressure

; set parameters for individual problem;  default ion=0, nocon=0
; ion can be set in eit_dem_tool, such as he2, si11, fe14, etc.
; nocon =  0 or 1 as desired
nocon=0

if n_elements(instr) eq 0 then instr='eit'
if not keyword_set(wmin) then wmin = 165
if not keyword_set(wmax) then wmax = 350
if n_elements(delta_wave) ne 0 then begin
    if not keyword_set(irradiance) then message,'Must set IRRADIANCE with '+$
            'DELTA_WAVE keyword'
    nbins = (wmax-wmin)/float(delta_wave)
    dn_s = fltarr(nbins)
    dum = dn_s
endif else dn_s = 0.
if n_elements(coefs) eq 0 then use_kc = 0 else use_kc = 1
bad = -1
elog = alog10(1./0.4343)
pixrad = (214.94*960./eit_pixsize())^2
if keyword_set(nocon) then begin
      cont = 0
    app = 'nocon'
endif else begin
    cont= 1
    app = 'con'
endelse
sfile = 'radiance_unit_dem_'+app+'.save'
case 1 of
	;file_test(concat_dir('$SSW_EIT_RESPONSE',sfile) : restore,concat_dir('$SSW_EIT_RESPONSE',sfile) 
	file_test(sfile): restore,sfile
	else: begin
		mk_eit_spec,cont=cont,he_fac=he_fac,abund=abund,eqion=eqion,pressure=pressure,file=sfile
		restore,sfile
	   endelse
endcase
if keyword_set(ion) then begin
	ion = ion+'_'
	sfile = ion + 'radiance_unit_dem_'+app+'.save'
	if not file_test(sfile) then mk_eit_spec,cont=cont,he_fac=he_fac,abund=abund,$
		eqion=eqion,pressure=pressure,/ion
	restore,sfile
	if ion eq 'he2_' then spec0 = spec0*he_fac
endif

case 1 of
;select only temperatures gt logT=4.6 for modeling EIT only DN values
   instr eq 'eit': begin
       instr_par = eit_parms(waves,wave,filter)
       ok = where_array(round(10*intemp),round(10*alog10(temp)))
       temp = temp(ok)
       spec0 = spec0(*,ok)
    end
;   instr eq 'lasco': instr_par = lasco_parms()
   else: instr_par = 1.
endcase

if n_elements(trange) eq 0 then trange = [0,n_elements(intemp)-1]
wrange = where(waves ge wmin and waves le wmax)
waves = waves(wrange)
spec0 = spec0(wrange,*)

sz=size(dem)
if sz(0) eq 1 then dem_map = dem + intemp + elog else begin
    dem_map = dem
    for i=0,sz(3)-1 do dem_map(0,0,i) = dem_map(*,*,i) + intemp(i) + elog
endelse

;compute line intensity for each temperature bin and sum
for i= trange(0),trange(1) do begin
    this_t = intemp(i)
    case 1 of
       use_kc    : begin	; used for Correction coefficient expansion
                   if this_t le x(1) then kcorr_t = coefs(*,*,0) * this_t + $
                             coefs(*,*,1) else kcorr_t = coefs(*,*,2)*this_t*$
                                     this_t+coefs(*,*,3)*this_t+coefs(*,*,4)

                   if this_t eq 6.3 then k63 = kcorr_t
                   if this_t eq 6.4 then kcorr_t = temporary(kcorr_t) > $
                                                      (k63 - 3) < (k63 + 0.6)
                   if this_t eq 6.5 then begin
                       bad = where(kcorr_t ge 2)
                       kcorr_t = temporary(kcorr_t) > (k63 - 6) < (k63 + 1.)
                   endif
                   dem_t = dem_map(i) + (kcorr_t < 4)
                   end
       sz(0) eq 1: dem_t = dem_map(i)	; input linear DEM function 
       else      : begin		; input DEM map, optional indices
                   dem_t = dem_map(*,*,i)
                    if n_elements(xyidx) ne 0 then begin
                       ndem = dem_t
                       dem_t = ndem*0
                       dem_t(xyidx(0,*),xyidx(1,*)) = ndem(xyidx(0,*),xyidx(1,*))
                    endif
                   end
    endcase

    if n_elements(delta_wave) eq 0 then begin
        dum = total(instr_par * spec0(*,i)) * 10^dem_t
        if bad(0) ne -1 then dum(bad) = 0.
        if keyword_set(irradiance) then dum = total(dum)/pixrad
    endif else begin
        for j = 0,nbins - 1 do dum(j) = total(spec0(j*delta_wave:j*delta_wave +$
                                  delta_wave-1,i))
        if bad(0) ne -1 then dem_t(bad) = 0.
        dum = (dum * total(10^dem_t))/pixrad
    endelse
; Simpson's rule integration
    if i eq trange(0) and (trange(1) - trange(0) gt 1) then dum = dum * 0.5
    if i eq trange(1) and (trange(1) - trange(0) gt 1) then dum = dum * 0.5
    dn_s = dn_s + dum*0.1
endfor
return,dn_s
end
