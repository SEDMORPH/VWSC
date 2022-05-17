;;** INPUT: 
; pcs_data supercolour amplitudes fltarr(namplitudes,ngal)
; pcerr    supercolour errors
; z        redshift (to check age of Universe)
; mag_abs  Absolute magnitude of galaxy at wave_norm
; pcs_model  Model supercolour amplitudes
; norm_model Normalisation of model
; wave_norm  Wavelength at which SEDs are normalised during PCA

;** OUTPUT
; lgm log stellar mass, current i.e. accounting for recyling


;** OPTIONAL OUTPUT
; sfr star formation rate, averaged over last 10^8 years (as Magphys does)
; lgssfr specific SFR, averaged over last 10^8 years
; age_wr light weighted age
; lgm0 log stellar mass, total i.e. integral SFR dt

FUNCTION VWSC_FITMASS, pcs_data, pcerr, z, norm_data, pcs_model, norm_model, params, wave_norm, lgm0=lgm0, lgsfr=lgsfr, lgssfr=lgssfr, lgage_wr=lgage_wr, f_burst=f_burst,tau_av=tau_av,z_met=z_met

  ;;-- basic parameters
  c_in_aa = 2.99792d18                  ;in [Ang/s]
  AB0 = 2.5*alog10(4*!DPI*3.0857^2*1d5/3.839)
  
  ngal = n_elements(z)
  nmod = n_elements(norm_model)

  ;;-- model  
  f_nu = norm_model*(double(wave_norm)^2/c_in_aa) ;in Lsol/Hz
  ABmag = AB0 - 2.5* alog10(f_nu) -48.6
  
  mag_sun = 4.55                ;This is at 1micron, guess between z and J on http://mips.as.arizona.edu/~cnaw/sun.html
                                ;doesn't matter as put in and take out again

  ;;-- data
  norm_nu_rest = norm_data*double(wave_norm)^2/c_in_aa ;=norm_nu_obs/(1+zphot)
  mag_app = janskytomag(norm_nu_rest)               ;apparent magnitude at 1micron (rest)  
  dist = lumdist(z)                                    ;luminosity distance in Mpc
  mag_abs = mag_app -5*alog10(dist)-25           ;absolute magnitude at 1micron (rest)


  ;;-- parameters to be fitted
  ;;params_fit = replicate({lgML:0d, lgML0:0d, sfr:0d, age_wr:0d},nmod)
  params_fit = replicate({lgML:0d, lgML0:0d, lgsfr:0d, lgssfr:0d, lgage_wr:0d, f_burst:0d, tau_av:0d, z_met:0d},nmod)
  params_fit.lgML = alog10(params.mstr1*10.^(0.4*(ABmag-mag_sun)))  ;current
  params_fit.lgML0 = alog10(params.mstr0*10.^(0.4*(ABmag-mag_sun))) ;total
;  params_fit.sfr = alog10(params.ftot[0]/1d6)
;  ;ftot[0] is fraction of total current stellar mass created in last
;  1d6 years 

  params_fit.lgSSFR = alog10(params.sfrav[2]/params.mstr1)          ;sfrav[2] is model SFR averaged over last 10^8 years: SFR for a 1Msol galaxy
  ;; SFR = SFRav/M_curr_model * M_curr_data. But we don't know
  ;; M_curr_data.
  ;; So SFR = SFRav/M_curr_model * (M/L)_curr_model / 10d^0.4(Mabs_data - Msol)
  params_fit.lgSFR = alog10((params.sfrav[2]/params.mstr1) *  10d^params_fit.lgML) ;sSFR of model * M/L_curr of model, then use Mabs below to get SFR

  params_fit.lgage_wr = alog10(params.age_wr)
;; added ---
  params_fit.f_burst = (params.fburst[3])
  params_fit.tau_av = (params.tauv0)
  params_fit.z_met = (params.zmet)

;; note that sfr is the sfr for a 1Msol galaxy, need to x by mass of galaxy
;; nbins_pdf = [120,120,120,120]
  nbins_pdf = [120,120,120,120,120,120,120,120]


  ;;-- output arrays
  bestchisq = fltarr(ngal)
  bestmodelid = lonarr(ngal)
  lgm = replicate({per16:0.0,per50:0.0,per84:0.0,bestmodelid:0l,bestchisq:0d},ngal)
;;  lgm0 =  (sfr = (age_wr = replicate({per16:0.0,per50:0.0,per84:0.0},ngal)))
  lgm0 =  (lgsfr = (lgssfr = (lgage_wr = (f_burst = (tau_av = (z_met = replicate({per16:0.0,per50:0.0,per84:0.0},ngal)))))))


  print, 'fitting massses: ', ngal
  for i=0L,ngal -1 do begin
     if total(pcs_data[0,i]) eq 0 then continue
     if i/1000 eq i/1000. then print,i
     age_universe = galage(z[i],1000,/silent)
     ind_allowed = where(params.tform lt age_universe)
     vwsc_pdfcalc, pcs_data[*,i], pcerr[*,i], pcs_model[*,ind_allowed], params_fit[ind_allowed], nbins_pdf, pdf, xx, $
                   stats_pdf=stats_pdf,stats_logL=stats_logL

     lgm[i].bestchisq = stats_logL.minchisq
     lgm[i].bestmodelid = ind_allowed[stats_logL.bestmodelid]

     lgML =  [stats_pdf[0].per16,stats_pdf[0].per50, stats_pdf[0].per84]
     lgm[i].per16 = lgML[0] - 0.4*(mag_abs[i]-mag_sun)
     lgm[i].per50 = lgML[1] - 0.4*(mag_abs[i]-mag_sun)
     lgm[i].per84 = lgML[2] - 0.4*(mag_abs[i]-mag_sun)

     lgML0 =  [stats_pdf[1].per16,stats_pdf[1].per50, stats_pdf[1].per84]
     lgm0[i].per16 = lgML0[0] - 0.4*(mag_abs[i]-mag_sun)
     lgm0[i].per50 = lgML0[1] - 0.4*(mag_abs[i]-mag_sun)
     lgm0[i].per84 = lgML0[2] - 0.4*(mag_abs[i]-mag_sun)

;     sfr[i].per16 = stats_pdf[2].per16
;     sfr[i].per50 = stats_pdf[2].per50
;     sfr[i].per84 = stats_pdf[2].per84

     lgssfr[i].per16 = stats_pdf[3].per16
     lgssfr[i].per50 = stats_pdf[3].per50
     lgssfr[i].per84 = stats_pdf[3].per84

     lgsfr[i].per16 = stats_pdf[2].per16 - 0.4*(mag_abs[i]-mag_sun)
     lgsfr[i].per50 = stats_pdf[2].per50 - 0.4*(mag_abs[i]-mag_sun)
     lgsfr[i].per84 = stats_pdf[2].per84 - 0.4*(mag_abs[i]-mag_sun)

     lgage_wr[i].per16 = stats_pdf[4].per16
     lgage_wr[i].per50 = stats_pdf[4].per50
     lgage_wr[i].per84 = stats_pdf[4].per84

;; added ---

     f_burst[i].per16 = stats_pdf[5].per16
     f_burst[i].per50 = stats_pdf[5].per50
     f_burst[i].per84 = stats_pdf[5].per84

     tau_av[i].per16 = stats_pdf[6].per16
     tau_av[i].per50 = stats_pdf[6].per50
     tau_av[i].per84 = stats_pdf[6].per84

     z_met[i].per16 = stats_pdf[7].per16
     z_met[i].per50 = stats_pdf[7].per50
     z_met[i].per84 = stats_pdf[7].per84

 
  endfor

  return, lgm

END
