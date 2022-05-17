PRO VWSC_PHOTOGRID, lambda, spec, minz, maxz, dz, wavemin, wavemax, wave_norm,filterset, dir_filters, gridfile


;;-- set up redshift array
  nredshift = round((maxz-minz)/dz)+1
  redshift = findgen(nredshift)*dz + minz 

;;-- filters
  readcol, filterset, filterfiles,ll_eff,filternames,form='(A, F,A)',/silent ;list of filter files
  nband = n_elements(ll_eff)
  
;;-- Read in model flux arrays
  ngal = (size(spec,/dim))[1]

  nbin = nband*nredshift
  flux_super = dblarr(nbin,ngal)


;;-- for each band convolve at all redshifts required
;;-- output from vw_conv_filter is in f_lambda_obs
;;-- flux_super is in f_lambda_rest
  nn=0L
  for j=0,nband-1 do begin
     
     readcol, dir_filters+filterfiles[j],lam_fil,fil,skipline=1,/silent
;;     if j le 11 then lam_fil = lam_fil*10000
     print, 'convolving band ',filternames[j]

     ind_j = indgen(nredshift)+nn 
     
     flux_super[ind_j,*] =  vwsc_conv_filter(lambda,spec,lam_fil,fil,redshift)
        
     ;; all fluxes are *rest frame* 
     flux_super[ind_j,*] = flux_super[ind_j,*]*rebin((1+redshift),nredshift,ngal) ;get rid of redshift effect for easier visualisation
     
     nn+=nredshift
             
  endfor
  
;; store supersampled wavelength, redshift and filtername arrays
  wave_rest_super = (wave_obs_super = (z_super = fltarr(nbin)))
  filternames_super = strarr(nbin)
  
  for k=0,nband-1 do begin
     wave_obs_super[k*nredshift:(k+1)*nredshift-1] = ll_eff[k]
     wave_rest_super[k*nredshift:(k+1)*nredshift-1] = ll_eff[k]/(1+redshift)
     filternames_super[k*nredshift:(k+1)*nredshift-1] = filternames[k]
     z_super[k*nredshift:(k+1)*nredshift-1] = redshift
  endfor

;; cut down to minimum and maximum wavelengths required
  ind_wave = where(wave_rest_super gt wavemin and wave_rest_super lt wavemax,npix)
  flux_super = flux_super[ind_wave,*]
  wave_rest_super = wave_rest_super[ind_wave]
  wave_obs_super = wave_obs_super[ind_wave]
  filternames_super = filternames_super[ind_wave]
  z_super = z_super[ind_wave]
  
;; calculate normalisation
  norm_super = fltarr(ngal)
  for i=0L,ngal-1 do begin
     ind2 = sort(wave_rest_super)     
     linterp, wave_rest_super[ind2], flux_super[ind2,i],wave_norm,tmp
     norm_super[i] = tmp
  endfor

  save, flux_super, norm_super, ind_wave, wave_obs_super, wave_rest_super, filternames_super, z_super, file = gridfile

END

