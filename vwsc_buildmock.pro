;;*** AIM: Build a mock dataset from the models, using a small number of all the available filters and random redshifts

PRO VWSC_BUILDMOCK, gridfile, evecfile, dir_figs, filterset,ind_filt, wave_norm, minz, maxz, dz, zphot_err, namplitudes, mockfile

;; this is the super sampled array, in bins of dz
  restore, gridfile             ;flux_super, norm, ind_wave, wave_obs_super, wave_rest_super, filternames_super, z_super, 
  npix = n_elements(wave_rest_super)
  nobj = (size(flux_super,/dim))[1]

  restore, evecfile

;;-- filters
  readcol, filterset, filterfiles,ll_eff,filternames,form='(A, F,A)',/silent ;list of filter files
  nband = n_elements(ll_eff)

;;------------------------------------------------------------------
;; Create mock  dataset, setting to zero filters we don't want
;; This represents an input data array i.e. nband x ngal
;;------------------------------------------------------------------
  flux_super_err = flux_super/10.                           
  seed = 2340975
  zmod = randomu(seed,nobj)*(maxz-minz)+minz   ;assign galaxies a redshift
  zmod_phot = zmod+randomn(seed,nobj)*zphot_err ;> minz < maxz
  
;; this is the mock observations to match catalog
  ind_ll_mock = ind_filt
  ll_eff_mock = ll_eff[ind_ll_mock]
  filtername_mock = filternames[ind_ll_mock]
  nband_mock = n_elements(ind_ll_mock)

  flux_mock = (flux_err_mock = fltarr(nband,nobj))
  flux_mock_zphot = (flux_err_mock_zphot = fltarr(nband,nobj))

  for i=0l,nobj-1 do begin
     tmp = abs(zmod[i]-z_super) ; should give all filters in utilised rest wavelength range
     ind_zz = where(tmp eq min(tmp))     

     tmp = abs(zmod_phot[i]-z_super)
     ind_zphot = where(tmp eq min(tmp))

     for j=0,nband_mock-1 do begin
        ind_ll = where(wave_obs_super eq ll_eff_mock[j])
        ind = vw_setintersection(ind_ll,ind_zz)
        if ind[0] ne -1 then begin
           flux_mock[ind_ll_mock[j],i] = flux_super[ind,i]
           flux_err_mock[ind_ll_mock[j],i] = flux_super_err[ind,i]     
        endif
        ind = vw_setintersection(ind_ll,ind_zphot)
        if ind[0] ne -1 then begin
           flux_mock_zphot[ind_ll_mock[j],i] = flux_super[ind[0],i] ;sometimes get >1 match where zphot lies directly between 2 z bins
           flux_err_mock_zphot[ind_ll_mock[j],i] = flux_super_err[ind[0],i]     
        endif
     endfor
  endfor


;;------------------------------------------------------------------
;; Place mock-UDS data array into the full matrix, creating a sparsely
;; sampled matrix. 
;;------------------------------------------------------------------

;;-- the sparsely sampled matrix
  flux_super_mock = (flux_super_mock_err = fltarr(npix,nobj))
  flux_super_mock_zphot = (flux_super_mock_err_zphot = fltarr(npix,nobj))
  
  for i=0l,nobj-1 do flux_super_mock[*,i] = vwsc_fillflux(flux_mock[*,i],zmod[i],minz,maxz,dz,ll_eff,ind_wave)
  for i=0l,nobj-1 do flux_super_mock_err[*,i] = vwsc_fillflux(flux_err_mock[*,i],zmod[i],minz,maxz,dz,ll_eff,ind_wave)

  for i=0l,nobj-1 do flux_super_mock_zphot[*,i] = vwsc_fillflux(flux_mock_zphot[*,i],zmod_phot[i],minz,maxz,dz,ll_eff,ind_wave)
  for i=0l,nobj-1 do flux_super_mock_err_zphot[*,i] = vwsc_fillflux(flux_err_mock_zphot[*,i],zmod_phot[i],minz,maxz,dz,ll_eff,ind_wave)
  
;;-- plot number observed bins per galaxy
  project = uintarr(npix,nobj)
  project[where(flux_super_mock ne 0)] = 1
  project2 = uintarr(npix,nobj)
  project2[where(flux_super_mock_zphot ne 0)] = 1

  ps1,dir_figs+'vwsc_buildmock_nfilters.ps'
  hist = histogram(total(project,1),binsize=1,min=0,max=9,loc=xx) ; xx are starting locations for bins
  plot, xx,hist,psym=10,xtitle='Number of filters',ytitle='Number of galaxies',xr=[0,10] ;psym=10 puts center of bins at xx
  hist = histogram(total(project2,1),binsize=1,min=0,max=9,loc=xx)                        ; xx are starting locations for bins
  plot, xx,hist,psym=10,xtitle='Number of filters',ytitle='Number of galaxies',xr=[0,10],title='with Zphot' ;psym=10 puts center of bins at xx
  ps2


;;------------------------------------------------------------------
;;-- normalise 
;;------------------------------------------------------------------
  flux_super_mock_norm = (flux_super_mock_err_norm = fltarr(npix,nobj))
  flux_super_mock_zphot_norm = (flux_super_mock_err_zphot_norm = fltarr(npix,nobj))
 
  ind2 = sort(wave_rest_super)  ;correct normalisation
  norm_gappy = fltarr(nobj)

  for i=0L,nobj-1 do begin

     ;; approximate normalisation, to stop PCA falling over

     ;; for normgappy
     ind = where(project[ind2,i] eq 1) 
     linterp, wave_rest_super[ind2[ind]], flux_super_mock[ind2[ind],i],wave_norm, tmp
     flux_super_mock_norm[*,i] = flux_super_mock[*,i]/tmp[0]
     flux_super_mock_err_norm[*,i] = flux_super_mock_err[*,i]/tmp[0]
     
     ;; for zphot + normgappy
     flux_super_mock_zphot_norm[*,i] = flux_super_mock_zphot[*,i]/tmp[0]
     flux_super_mock_err_zphot_norm[*,i] = flux_super_mock_err_zphot[*,i]/tmp[0]
     
     norm_gappy[i] = tmp[0]     ;stored for plotting purposes
  endfor

;;------------------------------------------------------------------
;;-- PCA
;;------------------------------------------------------------------
 
  pcs_mock = vwpca_normgappy(flux_super_mock_norm,flux_super_mock_err_norm, evecs[*,0:namplitudes-1], meanarr,norm=norm_mock)
  norm_mock = norm_gappy*norm_mock ;include approximate normalisation done above
  pcs_mock_zphot = vwpca_normgappy(flux_super_mock_zphot_norm,flux_super_mock_err_zphot_norm, evecs[*,0:namplitudes-1], meanarr,norm=norm_mock_zphot)
  norm_mock_zphot = norm_gappy*norm_mock_zphot

  save, flux_mock, flux_mock_zphot, flux_err_mock, flux_err_mock_zphot, flux_super_mock_zphot, flux_super_mock_err_zphot, flux_super_mock, flux_super_mock_err, zmod, zmod_phot, pcs_mock, pcs_mock_zphot, norm_mock,norm_mock_zphot, file = mockfile


END
