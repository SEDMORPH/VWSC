PRO vwsc_catalog_cosmos2015, catfile, fluxfile, ind_filt, filterset, minz, maxz, dir_notes,klim,maskfile=maskfile,zfile=zfile


  splog, file=dir_notes+'vwsc_catalog_cosmos.splog'

;;------------------------------------------------------------------ 
;;-- input files 
;;------------------------------------------------------------------

;;-- filterfiles: needed for E(B-V) correction
  readcol, filterset, filterfiles,ll_eff,name,form='(A,F,A,X)',/silent ;list of filter files
  nband = n_elements(ll_eff)

  ind_select = where(strmatch(name, 'Kvis') eq 1) ;;Kvis for COSMOS, Kuk for UDS
  ind_irac = where(strmatch(name, '3p6') eq 1 or strmatch(name, '4p5') eq 1) ;need to add other filters if ever use them
  splog, 'Selection filter:', name[ind_select]
  splog, 'Irac filters: ', name[ind_irac]

;;-- COSMOS data file
  data = mrdfits(catfile,1,header)
  splog, '# galaxies in data file:',n_elements(data)

  zspec = data.zspec_zCOSMOSDR3
  zused = data.PHOTOZ ;start with photometric redshifts
  ind = where(zspec gt minz and zspec lt maxz, compl=compl)
  zused[ind] = zspec[ind] ;replace zCOSMOS redshifts where possible
  
  zspec_chandra = data.z_best_chandra
  ind = where(zspec_chandra gt minz and zspec_chandra lt maxz, compl=compl)
  zused[ind] = zspec_chandra[ind] ;replace AGN redshifts where possible (some may include photometric instead of spec)

  ind = where(zused ge minz and zused le maxz,ngal)
  data = data[ind]
  ind = where(data.TYPE eq 0 or data.TYPE eq 2,ngal) ;selects galaxies and X-ray sources (AGN)
  data = data[ind]
  splog, '# of galaxies in z range: ',ngal

;;-- get matching survey area ---------
;;  if n_elements(maskfile) ne 0 then begin
;;     mask = mrdfits(maskfile,0,hdr)
;;     area_pix = abs(sxpar(hdr,'CD1_1')*sxpar(hdr,'CD2_2'))
;;     area_survey = n_elements(where(mask eq 0))*area_pix
;;     splog, 'area survey (degrees): ',area_survey
;; endif

  area_survey = 1.38
  splog, 'area survey (degrees): ',area_survey


;;------------------------------------------------------------------
;;-- calculate fluxes
;;------------------------------------------------------------------

;; dust-corrected mags_total = mags_3" - dust_corr + mag_offset -
;;                                       systematic offset
;; dust_corr = EBV*dust extinction factor


dustfactor = [4.02,3.117,2.66,1.991,1.461,1.211,0.871,0.563,0.364,0.162,0.111] ; 27/06/2018 - EBV*dust extinction factor

  mag = fltarr(nband,ngal)
  mag[ind_filt[0],*] = data.B_MAG_APER3+data.offset-0.146-data.ebv*dustfactor[0]
  mag[ind_filt[1],*] = data.V_MAG_APER3+data.offset+0.117-data.ebv*dustfactor[1]
  mag[ind_filt[2],*] = data.R_MAG_APER3+data.offset+0.012-data.ebv*dustfactor[2]
  mag[ind_filt[3],*] = data.ip_MAG_APER3+data.offset-0.02-data.ebv*dustfactor[3]
  mag[ind_filt[4],*] = data.zpp_MAG_APER3+data.offset+0.084-data.ebv*dustfactor[4]
  mag[ind_filt[5],*] = data.Y_MAG_APER3+data.offset-0.001-data.ebv*dustfactor[5]
  mag[ind_filt[6],*] = data.J_MAG_APER3+data.offset-0.017-data.ebv*dustfactor[6]
  mag[ind_filt[7],*] = data.H_MAG_APER3+data.offset-0.055-data.ebv*dustfactor[7]
  mag[ind_filt[8],*] = data.Ks_MAG_APER3+data.offset+0.001-data.ebv*dustfactor[8]
  mag[ind_filt[9],*] = data.SPLASH_1_MAG+0.025-data.ebv*dustfactor[9]
  mag[ind_filt[10],*] = data.SPLASH_2_MAG+0.005-data.ebv*dustfactor[10]; irac mags are already total, so no aperture correction needed


  magerr = fltarr(nband,ngal)
  magerr[ind_filt[0],*] = data.B_MAGERR_APER3
  magerr[ind_filt[1],*] = data.V_MAGERR_APER3
  magerr[ind_filt[2],*] = data.R_MAGERR_APER3
  magerr[ind_filt[3],*] = data.ip_MAGERR_APER3
  magerr[ind_filt[4],*] = data.zpp_MAGERR_APER3
  magerr[ind_filt[5],*] = data.Y_MAGERR_APER3
  magerr[ind_filt[6],*] = data.J_MAGERR_APER3
  magerr[ind_filt[7],*] = data.H_MAGERR_APER3
  magerr[ind_filt[8],*] = data.Ks_MAGERR_APER3
  magerr[ind_filt[9],*] = data.SPLASH_1_MAGERR
  magerr[ind_filt[10],*] = data.SPLASH_2_MAGERR



;;-- apply magnitude cut
  ind = where(mag[ind_select,*] le klim and mag[ind_select,*] gt 0,ngal)
;  ind = where(mag[ind_select,*] le klim,ngal)
  splog, '# gal in K limit: ',klim, ngal
  data = data[ind]
  mag = mag[*,ind]
  magerr = magerr[*,ind]

;;-- convert to fluxes
  flux = magtojansky(mag)
  flux_err = flux-magtojansky(mag+magerr)


;;-- error floors: 5% for all and 20% for IRAC
  flux_err[ind_filt,*] >= 0.05*flux[ind_filt,*]    ;5% floor on all
  flux_err[ind_irac,*] >= 0.2*flux[ind_irac,*] ;irac bands: this is about right c.f. model residuals (figs_testing.pro)



;;------------------------------------------------------------------
;;-- output for PCA 
;;------------------------------------------------------------------


;;-- bins that get no weight in pca fit
  ID = data.NUMBER
  zspec = data.zspec_zCOSMOSDR3
  ra = data.ALPHA_J2000
  dec = data.DELTA_J2000
  zphot = data.PHOTOZ
  
  zphot_chandra = data.z_phot_chandra
  ind = where(zphot_chandra gt 0)
  zphot[ind]=zphot_chandra[ind] ;photometric redshifts for AGN for type=0 objects
  
  ind = where(zspec gt minz and zspec lt maxz, compl=compl)
  zused = zphot
  zused[ind] = zspec[ind]
  
  zspec_chandra = data.z_spec_compile
  ind = where(zspec_chandra gt minz and zspec_chandra lt maxz, compl=compl)
  zspec[ind] = zspec_chandra[ind]
  
  Kmag=mag[ind_select,*]
  
;;-- bins that get no weight in pca fit
  ind_bad = where(flux lt 0 or mag gt 99 or mag eq 0)
  flux[ind_bad] = 0.0
  flux_err[ind_bad]=0.0         

  
  save, ID, flux, flux_err, zphot, zspec, zused, mass, Kmag, ra, dec, area_survey, file = fluxfile

  splog,/close


END
