;+
; NAME:
;	VWSC_ABMAG_Z0
;
; PURPOSE:
;	 calculate and output to file rest-frame u, V, K, J and 1mum
;	 absolute magnitudes for stochastic burst library SEDs
;
;
; CALLING SEQUENCE:
;	vwsc_abmag_z0, lambda, sed, [filterfile], [ll_eff]
;
; INPUTS:
;	lambda: model wavelength (float array[Npix])
;	sed: model SEDs in Lsol/AA (float array[Npix,Nmod])
;
; OPTIONAL INPUT: 
;       filterfile: String containing full path to filterfile, in two-column format: lambda_AA flux  
;       ll_eff: Effective wavelength of filter
; 
; KEYWORD PARAMETERS:
;	WAVELENGTH: If no filterfile is provided, calculate AB mag at this fixed wavelength
;
;
; EXAMPLE:
;       VWSC_ABMAG_Z0, lambda, sed, filterfile, ll_eff
;       or
;       VWSC_ABMAG_Z0, lambda, sed, wavelength = X
;
; MODIFICATION HISTORY:
; 	Written by:	Vivienne Wild 2013
;	May, 2014	Included into VWSC package
;-


FUNCTION VWSC_ABMAG_Z0, lambda, sed, filterfile, ll_eff, wavelength=wavelength

  @vwsc_constants.inc
  
  if n_elements(filterfile) ne 0 and keyword_set(wavelength) then message, 'Provide filterfile OR wavelength'

  ;; see USEFUL/ASTRO/ABmags.txt
  ;; f_nu needs converting into a flux density [erg/s/cm^2/Hz] for object at 10pc
  ;; 10pc = 3.0857e19cm; Lsol = 3.839e33erg/s
  ;; flux density = F_nu / (4 * !PI * d^2)
  AB0 = 2.5*alog10(4*!DPI*3.0857^2*1d5/3.839)

  if n_elements(filterfile) ne 0 then begin
     ;;-- filter profile
     readcol, filterfile,lam_fil,fil,/silent
     
     ;;-- magnitude
     flux = vwsc_conv_filter(lambda,sed,lam_fil,fil,0.0)
     f_nu = flux*(ll_eff^2/c_in_AA)[0] ;in Lsol/Hz
     ABmag = AB0 - 2.5* alog10(f_nu) -48.6
     
  endif else if keyword_set(wavelength) then begin

     ;;-- fixed wavelength
     ind = (where(lambda gt wavelength))[0]
     flux = sed[ind,*] 
     f_nu = flux*(wavelength^2/c_in_AA)  ;in Lsol/Hz
     ABmag = reform(AB0 - 2.5* alog10(f_nu) -48.6)
  endif else message, 'Please provide a filterfile, or specify a wavelength'

return, ABmag

END
