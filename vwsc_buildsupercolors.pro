FUNCTION FUNC_PHOTOGRID_READFILES, filename, nfiles, ngal=ngal, params=params

  
  BC03_READ_STOCH, filename, params, spec, wave,nfile=nfiles
  
  npix = n_elements(wave)
  ngal = (size(spec,/dim))[1]
  struct = {sed:fltarr(npix,ngal), lambda:fltarr(npix)}
  
  struct.sed = spec
  struct.lambda = wave
  
  return, struct
END



PRO VWSC_BUILDSUPERCOLORS, outfile, filterset, bc03fileid, minz, maxz, dz, wave_norm, wavemin_rest, wavemax_rest, namplitudes,zphot_err,nstoc=nstoc

  dir_filters = !dataDIR+'FILTERS/'


;;-- set up redshift array
  splog,'minz maxz  dz: ', minz, maxz, dz
  nredshift = round((maxz-minz)/dz)+1
  redshift = findgen(nredshift)*dz + minz 

;;-- filters
  readcol, filterset, filterfiles,ll_eff,filternames,form='(A, F,A)',/silent ;list of filter files
  nband = n_elements(ll_eff)
  
;;-- Read in model flux arrays
  if not(keyword_set(nstoc)) then nstoc=4 ;between 1 and 4 depending on how many models you want to use
  sed = func_photogrid_readfiles(bc03fileid, nstoc, ngal=ngal ,params=params) ;ngal is output here

  nbin = nband*nredshift
  flux_super = dblarr(nbin,ngal)

  splog, 'number of seds in file:',ngal

;;-- Luminosity distance for apparent magnitudes
;  Dl = lumdist(redshift,/silent)             ;Mpc
  c=2.99792d18                               ;in [Ang/s]

;;-- read in filter file - 300 is sufficient
  filter = replicate({lam_fil:fltarr(300), fil:fltarr(300),n:0L},nband)
  for j=0,nband -1 do begin
     readcol, dir_filters+filterfiles[j],lam_fil,fil,/silent
     n = n_elements(fil)
     filter[j].lam_fil[0:n-1] = lam_fil
     filter[j].fil[0:n-1] = fil
     filter[j].n = n
  endfor


;;-- for each band convolve at all redshifts required
;;-- output from vw_conv_filter is in f_lambda_obs
;;-- flux_super is in f_lambda_rest
  nn=0L
  for j=0,nband-1 do begin
     
     print, 'convolving band ',filternames[j]

     ind_j = indgen(nredshift)+nn ;[nn:nn+nredshift-1]
     
     flux_super[ind_j,*] =  vwsc_conv_filter(sed.lambda,sed.sed,filter[j].lam_fil[0:filter[j].n-1],filter[j].fil[0:filter[j].n-1],redshift)
        
     ;; all fluxes are *rest frame* 
     flux_super[ind_j,*] = flux_super[ind_j,*]*rebin((1+redshift),nredshift,ngal) ;get rid of redshift effect for easier visualisation
     
     nn+=nredshift
             
  endfor
  
  

  nbin = (size(flux_super,/dim))[0]
  nobj = (size(flux_super,/dim))[1]

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
  ind_wave = where(wave_rest_super gt wavemin_rest and wave_rest_super lt wavemax_rest,npix)
  flux_super = flux_super[ind_wave,*]
  wave_rest_super = wave_rest_super[ind_wave]
  wave_obs_super = wave_obs_super[ind_wave]
  filternames_super = filternames_super[ind_wave]
  z_super = z_super[ind_wave]


;;-- normalise
  norm_models  = fltarr(ngal)
  flux_super_norm = flux_super
  for i=0,ngal-1 do begin
     ind2 = sort(wave_rest_super) ;correct normalisation
     linterp, wave_rest_super[ind2], flux_super[ind2,i],wave_norm,tmp
     norm_models[i] = tmp
     flux_super_norm[*,i] = flux_super[*,i]/norm_models[i]
  endfor


;;------------------------------------------------------------------
;;-- PCA
;;------------------------------------------------------------------
  evecs = vwpca(flux_super_norm,namplitudes,variance,pcs,meanarr,evals,/svd)
  splog, 'variances: ',variance[0:namplitudes-1]
  splog, 'total variance in evecs: ',total(variance[0:namplitudes-1])

;;------------------------------------------------------------------
;; Create mock dataset
;;------------------------------------------------------------------
  flux_super_err = flux_super/10.                           
  seed = 2340975
  zmod = randomu(seed,nobj)*(maxz-minz)+minz   ;assign galaxies a redshift
  if zphot_err ne 0.0 then z_normgappy = zmod+randomn(seed,nobj)*zphot_err else z_normgappy = zmod

  flux_uds = (flux_err_uds = fltarr(nband,nobj))

  for i=0l,nobj-1 do begin
     tmp = abs(z_normgappy[i]-z_super) ; should give all filters in utilised rest wavelength range
     ind_zz = where(tmp eq min(tmp))     

     for j=0,nband-1 do begin
        ind_ll = where(wave_obs_super eq ll_eff[j])
        ind = setintersection(ind_ll,ind_zz)
        if ind[0] ne -1 then begin
           flux_uds[j,i] = flux_super[ind,i]
           flux_err_uds[j,i] = flux_super_err[ind,i]     
        endif
     endfor
  endfor

  ;;-- the sparsely sampled matrix
  flux_super_uds = (flux_super_uds_err = fltarr(npix,nobj))
  
  for i=0l,nobj-1 do flux_super_uds[*,i] = vwsc_fillflux(flux_uds[*,i],z_normgappy[i],redshift,ll_eff,ind_wave)
  for i=0l,nobj-1 do flux_super_uds_err[*,i] = vwsc_fillflux(flux_err_uds[*,i],z_normgappy[i],redshift,ll_eff,ind_wave)

  project = uintarr(npix,nobj)
  project[where(flux_super_uds ne 0)] = 1

  flux_super_uds_norm = (flux_super_uds_err_norm = fltarr(npix,nobj))
  norm_gappy = fltarr(nobj)
  for i=0,nobj-1 do begin
     ind = where(project[ind2,i] eq 1) ;approximate normalisation
     linterp, wave_rest_super[ind2[ind]], flux_super_uds[ind2[ind],i],wave_norm, tmp
     norm_gappy[i] = tmp

     ;; for normgappy
     flux_super_uds_norm[*,i] = flux_super_uds[*,i]/norm_gappy[i]
     flux_super_uds_err_norm[*,i] = flux_super_uds_err[*,i]/norm_gappy[i]
  endfor
 
  pcs_normgappy = vwpca_normgappy(flux_super_uds_norm,flux_super_uds_err_norm, evecs[*,0:namplitudes-1], meanarr)
  
;;------------------------------------------------------------------
;;-- output
;;------------------------------------------------------------------
;;-- needed for udsz test
  evecs[*,0] = -evecs[*,0] 
  evecs[*,1] = -evecs[*,1]
  evecs[*,2] = -evecs[*,2]
  pcs[0,*] = - pcs[0,*]
  pcs[1,*] = - pcs[1,*]
  pcs[2,*] = - pcs[2,*]
  pcs_normgappy[0,*] = - pcs_normgappy[0,*]
  pcs_normgappy[1,*] = - pcs_normgappy[1,*]
  pcs_normgappy[2,*] = - pcs_normgappy[2,*]

  save, params, evecs, meanarr, pcs,  variance, evals, norm_models, pcs_normgappy, z_normgappy, flux_super, redshift, ll_eff,ind_wave, wave_obs_super, wave_rest_super, filternames_super, z_super, file = outfile



END
