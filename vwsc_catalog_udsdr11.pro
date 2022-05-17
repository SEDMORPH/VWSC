PRO vwsc_catalog_udsdr11, catfile, fluxfile, ind_filt, filterset, minz, maxz, dir_notes,klim,maskfile=maskfile,area_survey=area_survey


  splog, file=dir_notes+'vwsc_catalog_udsdr11.splog'

;;------------------------------------------------------------------ 
;;-- input files 
;;------------------------------------------------------------------

;;-- filterfiles: needed for E(B-V) correction
  readcol, filterset, filterfiles,ll_eff,form='(A,F,X,X)',/silent ;list of filter files
  nband = n_elements(ll_eff)

;;-- UDS data file
  data = mrdfits(catfile,1,header)
  splog, '# galaxies in data file 1:',n_elements(data)

;;-- remove junk
  ;; these extra galaxies are PSBs and have good spectra, but are in
  ;; Y-band mask, so want to get their SCs. ALSO SET YBAND FLUX = 0.0
  ;; BELOW! 
  ;; ind = where((data.Best_galaxies eq 'T' or data.ID_1 eq '152227' or data.ID_1 eq '153502' or data.ID_1 eq '162358') and zused ge minz and zused le maxz,ngal)
  ind_best = where(data.Best_galaxies eq 'T')
  ind_noY = where(data.Best_galaxies eq 'F' and data.XTALK_DR11 eq 'F' and data.stars_final eq 'F' and (data.maskflag eq 8 or data.maskflag eq 9))

  ind = [ind_best, ind_noY]
  ngal = n_elements(ind)
  
  data = data[ind]

  splog, '# of Best galaxies: ',n_elements(ind_best)
  splog, 'Extra with bad Y-band: ',n_elements(ind_noY)
  splog, 'total number of galaxies:', ngal
  
;;------------------------------------------------------------------
;;-- remove galaxies outside redshift range
;;------------------------------------------------------------------

;use zspec in place of zphot for trimming data...

  zused = data.z_m2             ; see email conversation with Aaron and Omar about z_p vs z_m2 in Oct 2019
                                ; z_m2 gives a bunch of outliers below the MS, but removes "plume" going
                                ; through the PSB region
  zspec = data.z_spec
  ind = where(zspec ge 0)
  zused[ind] = zspec[ind]

  ind_z = where(zused ge minz and zused le maxz,ngal)
  data = data[ind_z]
  
  zused = -1
  ind_best = -1
  ind_noY_zrange = vw_setintersection(ind_noY,ind_z)
  ind_noY = -1
  splog, '# of objects in z range: ',ngal


;;-- get matching survey area

  if n_elements(maskfile) ne 0 then begin
     mask = mrdfits(maskfile,0,hdr)
     area_pix = abs(sxpar(hdr,'CD1_1')*sxpar(hdr,'CD2_2'))
     area_survey = n_elements(where(mask eq 0))*area_pix
     splog, 'area survey (degrees): ',area_survey
     splog, '# of galaxies in mask: ',ngal
  endif


;;------------------------------------------------------------------
;;-- calculate fluxes
;;------------------------------------------------------------------


  mag = fltarr(nband,ngal)
  mag[ind_filt[0],*] = data.BMAG_20 ;dust and aperture correction already applied
  mag[ind_filt[1],*] = data.VMAG_20
  mag[ind_filt[2],*] = data.RMAG_20
  mag[ind_filt[3],*] = data.iMAG_20
  mag[ind_filt[4],*] = data.zMAG_20
  mag[ind_filt[6],*] = data.JMAG_20
  mag[ind_filt[7],*] = data.HMAG_20
  mag[ind_filt[8],*] = data.KMAG_20
  mag[ind_filt[5],*] = data.YMAG_20
  mag[ind_filt[9],*] = data.I1MAG_20
  mag[ind_filt[10],*] = data.I2MAG_20

  magerr = fltarr(nband,ngal)
  magerr[ind_filt[0],*] = data.BMAGERR_20
  magerr[ind_filt[1],*] = data.VMAGERR_20
  magerr[ind_filt[2],*] = data.RMAGERR_20
  magerr[ind_filt[3],*] = data.iMAGERR_20
  magerr[ind_filt[4],*] = data.zMAGERR_20
  magerr[ind_filt[6],*] = data.JMAGERR_20
  magerr[ind_filt[7],*] = data.HMAGERR_20
  magerr[ind_filt[8],*] = data.KMAGERR_20
  magerr[ind_filt[5],*] = data.YMAGERR_20
  magerr[ind_filt[9],*] = data.I1MAGERR_20 
  magerr[ind_filt[10],*] = data.I2MAGERR_20


;;-- apply magnitude cut
  ind = where(mag[ind_filt[8],*] le klim,ngal)
  splog, '# gal in K limit: ',klim, ngal
  data = data[ind]
  mag = mag[*,ind]
  magerr = magerr[*,ind]
  ind_noY_zrange_maglim = vw_setintersection(ind_noY_zrange,ind)
  
  
;;-- convert to fluxes
  flux = magtojansky(mag)
  flux_err =  flux-magtojansky(mag+magerr) ;approx!

;;-- error floors: 5% for all and 20% for IRAC
  flux_err[ind_filt,*] >= 0.05*flux[ind_filt,*]    ;5% floor on all
  flux_err[ind_filt[9:10],*] >= 0.2*flux[ind_filt[9:10],*] ;irac bands: this is about right c.f. model residuals (figs_testing.pro)


;;------------------------------------------------------------------
;;-- output for PCA 
;;------------------------------------------------------------------

  ID = data.ID
  zspec = data.z_spec
  ra = data.ALPHA_J2000             ; From Nottingham mother catalogue
  dec = data.DELTA_J2000
  zphot = data.z_m2
  ind = where(zspec ge minz and zspec le maxz, compl=compl)
  zused = zphot
  zused[ind] = zspec[ind]
  Kmag=data.KMAG_20

;;-- bins that get no weight in pca fit
  ind_bad = where(flux lt 0 or mag gt 99 or mag eq 0) ;  or flux/flux_err lt 6)
  flux[ind_bad] = 0.0
  flux_err[ind_bad]=0.0

;; mask Y-band in the "extra" galaxies
  if ind_noY_zrange_maglim[0] ne -1 then begin
     splog, 'setting Y-band to zero in bad Y-band regions', n_elements(ind_noY_zrange_maglim)
     flux[ind_filt[5],ind_noY_zrange_maglim] = 0.0
     flux_err[ind_filt[5],ind_noY_zrange_maglim]=0.0
  endif
  
  save, ID, flux, flux_err, zphot, zspec, zused, Kmag, ra, dec, area_survey, file = fluxfile

  splog,/close


END
