;;-- EXAMPLE WRAPPER FILE FOR BUILDING SUPERCOLOURS

PRO VWSC_WRAPPER_COSMOS2015, NOMAG=NOMAG, NOCAT=NOCAT, NOGRID=NOGRID, NOPCA=NOPCA, NOMOCK=NOMOCK,NOMASS=NOMASS,NOMODELPCS=NOMODELPCS
  
;;****************************************************************************************
;;-- EDIT VARIABLES BELOW UNTIL YOU REACH "RUNNING" SECTION
;;****************************************************************************************

;;------------------------------------------------------------------
;;-- Variables: if any of these are edited, need to change dir_data
;;   below and rebuild eigenbasis
;;------------------------------------------------------------------
  modeltag_pca = 'rand'         ; ID of stochastic burst library file used to create PCA components
  nfile_pca = 4                 ; From 1-4, the number of stochastic burst library files to create PCA components
  
  modeltag_fit = 'all'          ; ID of stochastic burst library file used to fit stellar masses, get vmax etc. e.g. 'rand' for rough fit 'all' for better stats
  nfile_fit = 4                 ; From 1-4, the number of stochastic burst library files to fit M*, get vmax etc. e.g. 1 for rough fit, 4 for better stats

  wave_norm = 10000             ; Wavelength at which the eigenvectors will be normalised

  minz = 0.5                    ; Minimum redshift 
  maxz = 2.5                   ; Maximum redshift
  dz = 0.01                     ; Redshift intervals for filter shifts
  zphot_err = 0.06              ; Approximate error on photo-zs, for making mock catalogue

  wavemin = 2500.               ; Minimum rest-frame wavelength for eigenvectors
  wavemax = 15000.              ; Maximum rest-frame wavelength for eigenvectors

  namplitudes = 3               ;Number of eigenvectors to calculate. Usually 3 is sufficient to capture most of the information for galaxy SEDs. 
  
;;------------------------------------------------------------------
;;-- Directories
;;------------------------------------------------------------------
  fieldID = 'cosmos2015_ultradeep_wavemin2500_27_06_2018'  ;;;;;;;;;;;;;;;;;;;;;;;;;;;; Unique data identifier for dataset analysed in this wrapper
  z_ID = '0p5z2p5_wavemin2500' ;;;;;;;unique eigenvector identifier (redshift range and wavelength range dependent)


  dataDIR = '/home/ppxakw/data/'
  dir = dataDIR+'UDSPSB/' 
  dir_data = dir+'EBASIS/'      ; EBASIS directory needs to be changed if any of the "variables" above or filterID are changed

  dir2 = dir+fieldID+'/'             ; to match fieldID. Figures and notes are stored separately for each data catalogue
  fieldID = '_'+fieldID

  dir_notes = dir2+'NOTES/'
  dir_figs = dir2+'FIGS/'
  dir_output = dir2+'OUTPUT/'   ;for public consumption

  if file_test(dir2) eq 0 then spawn, 'mkdir '+dir2
  if file_test(dir_notes) eq 0 then spawn, 'mkdir '+dir_notes
  if file_test(dir_figs) eq 0 then spawn, 'mkdir '+dir_figs
  if file_test(dir_output) eq 0 then spawn, 'mkdir '+dir_output

  dir_model = dataDIR+'bc03/bc03_lib/'

  dir_cats = dataDIR+'CATS/COSMOS/'

  if modeltag_fit eq 'all' then readsfrav = 1 else readsfrav=0 ;instruct read_stoch to read SFRav as part of parameter file

;;------------------------------------------------------------------
;;-- FILTERS
;;------------------------------------------------------------------

  dir_filters = dataDIR+'FILTERS/'
  filterID = '_cosmos2015'                           ; Filterset identifier
  filterset = dir_filters+'vwsc'+filterID+'.lis' ; Must match filename of file containing filter list

;;------------------------------------------------------------------
;;-- Output Files (only output if not already made)
;;------------------------------------------------------------------

  fluxfile = dir_data+'flux'+fieldID+filterID+'_'+z_ID+'.sav'                ;data
  gridfile = dir_data+'modelgrid_'+modelTag_pca+filterID+'_'+z_ID+'.sav'     ;grid for PCA models
  gridfile_all = dir_data+'modelgrid_'+modelTag_fit+filterID+'_'+z_ID+'.sav' ;grid for parameter fitting models
  superfluxfile = dir_data+'superflux'+fieldID+filterID+'_'+z_ID+'.sav'      ;data
  evecfile = dir_data+'evecs_'+modelTag_pca+filterID+'_'+z_ID+'.sav'         ;models

  mockfile = dir_data+'modelmock_'+modelTag_pca+filterID+fieldID+'.sav'         ;models, to match data
  massfile = dir_data+'mass_'+modelTag_pca+'_'+modelTag_fit+filterID+fieldID+'.sav' ;data, depends on models
  massfile2 = dir_data+'mass2_'+modelTag_pca+'_'+modelTag_fit+filterID+fieldID+'.sav' ;data, depends on models
  magfile = dir_data+'ABmag_z0_'+modelTag_pca+'_'+modelTag_fit+filterID+fieldID+'.sav' ;data, depends on models
  modelpcfile = dir_data+'modelpcs_'+modelTag_pca+'_'+modelTag_fit+filterID+'_'+z_ID+'.sav'  ;all models, evecs depend on modelTag
  vmaxfile = dir_data+'vmax_'+modelTag_pca+'_'+modelTag_fit+filterID+fieldID+'_'+z_ID+'.sav'

  outfile = dir_output + 'COSMOS2015_'+string(minz,form='(F0.1)')+'_'+string(maxz,form='(F0.1)')+fieldID+'.fits'
  outfile2 = dir_output + 'COSMOS2015_'+string(minz,form='(F0.1)')+'_'+string(maxz,form='(F0.1)')+fieldID+'_topcat_v2.fits'  
;;------------------------------------------------------------------
;;-- Data
;;------------------------------------------------------------------
; Note that the catalogue file needs to be built using purpose written code. See e.g. vwsc_catalog_eg.pro
; suggest the "eg" includes field and filterset information
  
  ind_filt = [0,1,2,3,4,6,7,8,9,13,14] ; Filter numbers available  in data file  B, V, R, i, z++, (yHSC), Y, j, h, k, (G, IA738, IA768,)  ch1, ch2
  ;filters in brackets not used in data projection


  ind_select = 9               ;index of selection filter (K)
  maxmag_survey = 24.1

  catfile = dir_cats+'COSMOS2015_UVista_cat_withzCOSMOS_withChandra.fits''


  if not(keyword_Set(NOCAT)) then vwsc_catalog_cosmos2015, catfile, fluxfile, ind_filt, filterset, minz, maxz, dir_notes,maxmag_survey,maskfile= dir_cats+'ultravista.mask.stars_regions.fits'

;;****************************************************************************************
;;-- RUNNING... 
;;****************************************************************************************
;;------------------------------------------------------------------
;;-- Build eigenvectors
;;------------------------------------------------------------------

;;-- Load model data files 
  VWSC_READ_STOCH, dir_model+modeltag_pca, params, spec, lambda,nfile=nfile_pca;,maxgal=100000,readsfrav=readsfrav

;;-- Build Supersampled Photometric Grid and save
  if not(keyword_set(NOGRID)) then VWSC_PHOTOGRID, lambda, spec, minz, maxz, dz, wavemin, wavemax, wave_norm, filterset, dir_filters,gridfile


;;-- PCA to obtain Eigenvectors of supersampled photometric grid
  if NOT(keyword_set(NOPCA)) then begin
     restore, gridfile
     npix = (size(flux_super,/dim))[0] & ngal = (size(flux_super,/dim))[1]
     evecs = VWPCA(flux_super/transpose(rebin(norm_super,ngal,npix)),namplitudes,variance,pcs,meanarr,evals,/svd)
     
     ;; ensure the top 3 eigenvectors always come out the same way up
     ;; this only works if inputting fairly similar models to those used in Wild+2014
     indmin = where(wave_rest_super eq min(wave_rest_super))
     indmax = where(wave_rest_super eq max(wave_rest_super))
     if evecs[indmin,0] lt evecs[indmax,0] then begin
        evecs[*,0] = -evecs[*,0]
        pcs[0,*] = -pcs[0,*]
     endif
     if evecs[indmin,1] gt evecs[indmax,1] then begin
        evecs[*,1] = -evecs[*,1]
        pcs[1,*] = -pcs[1,*]
     endif
     if evecs[indmin,2] lt evecs[indmax,2] then begin
        evecs[*,2] = -evecs[*,2]
        pcs[2,*] = -pcs[2,*]
     endif

     save, evecs, variance, pcs, meanarr, evals, wave_rest_super, file=evecfile
  endif


;;-- Build a mock dataset, with random redshifts and only the filters that at available in the real catalog
;;-- calculate SCs for this mock dataset
  if NOT(keyword_set(NOMOCK)) then VWSC_BUILDMOCK, gridfile, evecfile, dir_figs, filterset,ind_filt, wave_norm,minz, maxz, dz, zphot_err, namplitudes, mockfile
  

;;-- Figures
;;if NOT(keyword_set(NOMOCK)) then VWSC_FIGS_MODELS, dir_figs, dir_notes, gridfile, mockfile, evecfile, filterset, lambda, spec, params,namplitudes, wave_norm
;;VWSC_FIGS_MODELS, dir_figs, dir_notes, gridfile, mockfile, evecfile, filterset, lambda, spec, params,namplitudes, wave_norm


;;-- apply to full set of models, for parameter fitting
  if not(keyword_set(NOMODELPCS)) then begin
     VWSC_READ_STOCH, dir_model+modeltag_fit, params_models, spec, lambda,nfile=nfile_fit,maxgal=100000,readsfrav=readsfrav
     VWSC_PHOTOGRID, lambda, spec, minz, maxz, dz, wavemin, wavemax, wave_norm, filterset, dir_filters,gridfile_all
     restore, gridfile_all
     restore, evecfile
     npix = (size(flux_super,/dim))[0] & nmod = (size(flux_super,/dim))[1]
     pcs_models =  (flux_super/transpose(rebin(norm_super,nmod,npix)) - rebin(meanarr,npix,nmod)) ## TRANSPOSE(evecs[*,0:namplitudes-1]) 
     save, pcs_models, params_models, file= modelpcfile
  endif


;;------------------------------------------------------------------
;;-- Apply eigenvectors to data 
;;------------------------------------------------------------------
  restore, fluxfile
  ngal = (size(flux,/dim))[1]
  
  pcs_data = (pcerr_data = fltarr(namplitudes,ngal))
  norm_data = fltarr(ngal)

  pcs_data = VWSC_PROJECTDATA(gridfile, evecfile, filterset, dir_figs,ID,flux, flux_err, zphot, namplitudes, minz, maxz, dz, chisq=chisq,ngood=ngood,pcerr=pcerr,norm=norm,superfluxfile=superfluxfile,maxplot=100)
  pcerr_data = pcerr
  norm_data = norm
  chisq_data = chisq
  ngood_data = ngood


  ind = where(zused ge minz and zused le maxz,compl=compl)
  pcs_data_zspec = (pcerr_zspec = fltarr(namplitudes,ngal))
  pcs_data_zspec[*,ind] = VWSC_PROJECTDATA(gridfile, evecfile, filterset,dir_figs,ID, flux[*,ind], flux_err[*,ind], zused[ind], namplitudes, minz, maxz, dz, chisq=chisq,ngood=ngood,pcerr=pcerr,norm=norm,maxplot=100)
  pcerr_zspec[*,ind] = pcerr
;;;;;  pcs_data_spec[*,compl]=-999.
  pcs_phot = pcs_data
  pcerr_phot = pcerr_data
  pcs_data[*,ind] = pcs_data_zspec[*,ind]
  pcerr_data[*,ind] = pcerr_zspec[*,ind]
;;;;  norm_data[*,ind] = norm
;;;;  chisq_data[*,ind] = chisq
;;;;  ngood_data[*,ind] = ngood

  ind = where(finite(pcs_data) eq 0) 
  if ind[0] ne -1 then begin
     pcs_data[ind]=-999
     pcerr_data[ind] = -999
  endif

;;------------------------------------------------------------------
;;-- Fit models to data using super-colours
;;------------------------------------------------------------------

;;-- calculate stellar masses from likelihood analysis of models
  if not(keyword_set(NOMASS)) then begin
     ;restore, evecfile
     restore, modelpcfile       ;pcs_models, params_models
     restore, gridfile_all      ;norm_super


     lgm = VWSC_FITMASS(pcs_data, pcerr_data, zused, norm_data, pcs_models, norm_super, params_models, wave_norm,lgsfr=lgsfr, lgssfr=lgssfr,age_wr=lgage_wr,f_burst=f_burst,tau_av=tau_av,z_met=z_met)

     mag_select_obs = reform(janskytomag(flux[ind_select,*])) < maxmag_survey ;allows for floating point errors going from flux to mag

    save, lgm, mag_select_obs,lgage_wr,lgsfr,lgssfr,f_burst,tau_av,z_met, file=massfile

  endif

  print, 'Fitting masses done'


;;-- Calculate standard AB magnitudes of models for comparison with traditional methods and for calculating stellar masses
; Williams et al. 2009 is the reference for UVJ colour diagram, they specify Bessell 1990 UV and Tokunaga 2002 mauna kea definition for J
; Bessell 1990 gives the filter functions in table2, but I downloaded them from here http://spiff.rit.edu/classes/phys440/lectures/filters/filters.html
; Effective wavelengths took A0V system table 3 bessell 1990

  if not(keyword_set(NOMAG)) then begin
     if file_test(massfile) eq 0 then begin
        print, 'Need to calculate bestmodel (massfile) to calculate model mags'
        stop
     endif
     restore, massfile          ;contains bestmodelid
     VWSC_READ_STOCH, dir_model+modeltag_fit, params_models, spec, lambda,nfile=nfile_fit,maxgal=100000,readsfrav=readsfrav
     npix =n_elements(lambda)
     ;; scale spectrum to median stellar mass

     spec_scaled=dblarr(npix,ngal)
     for i=0L,ngal-1 do spec_scaled[*,i]=spec[*,lgm[i].bestmodelid]*10d^lgm[i].per50/params_models[lgm[i].bestmodelid].mstr1


     ABmag_u_z0 = VWSC_ABMAG_Z0(lambda,spec_scaled,dir_filters+'Bessell/bess-u.pass',3659d)
     ABmag_V_z0 = VWSC_ABMAG_Z0(lambda,spec_scaled,dir_filters+'Bessell/bess-v.pass',5448d)
     ABmag_J_z0 = VWSC_ABMAG_Z0(lambda,spec_scaled,dir_filters+'UKIRT/wfcam_J.txt',12483.27d)
     ABmag_K_z0 = VWSC_ABMAG_Z0(lambda,spec_scaled,dir_filters+'UKIRT/wfcam_K.txt',22009.93d)
     ABmag_1mum_z0 = VWSC_ABMAG_Z0(lambda,spec_scaled,wavelength=1e4)
  
     save, ABmag_u_z0, ABmag_V_z0, ABmag_K_z0, ABmag_J_z0, ABmag_1mum_z0, file = magfile
  endif



;;-- the cuts have to be altered by hand :-( 

  lines = {red_sl:0.82, red_int:16.39, psb_sl:0.38,psb_int:11.37, psb2_sl:-0.26,psb2_int:3.94,dusty_sl:-0.18, dusty_int:-12.12,lomet_pc3lo:8, lomet_sl:0, lomet_pc1hi:10, green_cut:0, dust_cut:-19.5 }

  nzbin=4                       ;number of redshift bins to split data into for plotting


;;------------------------------------------------------------------
;;-- OUTPUT
;;------------------------------------------------------------------

  if file_test(massfile) ne 0 then restore, massfile
  if file_test(magfile) ne 0 then restore, magfile
  if file_test(magfile) ne 0 then uvj = transpose([[ABmag_u_z0],[ABmag_V_z0],[ABmag_J_z0]])


  VWSC_OUTPUT, pcs_data, pcerr_data, pcs_phot, pcerr_phot, ngood_data, norm_data, chisq_data, lines, fluxfile,outfile,outfile2,lgm=lgm,mass_scaled=lgm_scaled,lgage_wr=lgage_wr,lgsfr=lgsfr,lgssfr=lgssfr,uvj=uvj,f_burst=f_burst,tau_av=tau_av,z_met=z_met


END
