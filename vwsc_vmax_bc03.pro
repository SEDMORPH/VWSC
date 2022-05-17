;;*** VMAX_BC03.pro
;; AIM: return vmax for a galaxy
;; 
;; NOTES: this uses BC03 stochastic burst models to do K-correction,
;; have to tell it which file string and already have identified the
;; best fitting model for each galaxy
;;------------------------------------------------------------------

FUNCTION comoving, z

H0 = 70.                        ;km/s/Mpc
Olam = 0.7
Omas = 0.3
Otot = Omas+Olam
c = 3E5                     ;km/s


Dl = lumdist(z,H0=H0,Omega_m=Omas, Lambda0=Olam,/silent) ;Mpc

Rodr = (c/H0)*((1-Otot)*(1+z)^2+Olam+Omas*(1+z)^3)^(-0.5)

dV = Rodr*4*!PI*(Dl/(1.+z))^2

;print, z,dv,rodr,dl

return,dV

END

;;------------------------------------------------------------------

FUNCTION VWSC_VMAX_BC03, z_gal, appmag_observed_gal, appmag_obs_err, $
                    minmag_survey, maxmag_survey, minz_survey, maxz_survey, $
                    filename,bestmodel_gal, lam_fil,fil, $
                    z_kcorr=z_kcorr, nbin_kcorr = nbin_kcorr, spec=spec, wave=wave,nfile=nfile

common func, absmag, appmag, kcorr_igal,zgrid

ngal = n_elements(z_gal)

;;-- redshift to which K-correction is performed
;;-- use 0.0 to k-correct in "standard" way
if not(keyword_set(z_kcorr)) then z_kcorr = 0.0                   ;half way through survey? 

;;-- Number of redshift bins to use for K-correction
if not(keyword_set(nbin_kcorr)) then nbin_kcorr = 100

;;-- calulate total comoving volume of survey
vsurvey = QROMB('comoving',minz_survey,maxz_survey)

;;-- DM at edges of survey
DM_min = 5.*alog10(lumdist(minz_survey,/silent)*100000) ;10pc
DM_max = 5.*alog10(lumdist(maxz_survey,/silent)*100000) ;10pc

;;-- distance modulus for all the input galaxies
DM_gal = 5.*alog10(lumdist(z_gal,/silent)*100000)
ngal = n_elements(z_gal)

;;-- read in model SEDs (slow, best to pass in to debug)
if not(keyword_set(spec)) then begin
   if keyword_set(nfile) eq 0 then nfile=4
   VWSC_READ_STOCH, filename, params, spec, wave,nfile=nfile,maxgal=100000
   splog, 'number of seds in model file:',(size(spec,/dim))[1]
endif

;;------------------------------------------------------------------
;;-- calculate intrinsic absolute magnitude of galaxy in the survey
;;   selection band
;;------------------------------------------------------------------
;; all the constants cancel
;  f_nu = flux*(ll_eff[j]^2/c)[0]     ;in Lsol/Hz
;  ABmag_V_z0 = AB0 - 2.5* alog10(f_nu) -48.6

tmpmag_intrinsic  = -2.5*alog10(vwsc_conv_filter(wave,spec[*,bestmodel_gal],lam_fil,fil,z_kcorr))       
tmpmag_observed = fltarr(ngal)
for i=0L,ngal-1 do tmpmag_observed[i]  = -2.5*alog10(vwsc_conv_filter(wave,spec[*,bestmodel_gal[i]],lam_fil,fil,z_gal[i])) 
kcorr = (tmpmag_intrinsic - tmpmag_observed)                                                        

absmag_intrinsic_gal = (appmag_observed_gal+kcorr)-DM_gal ;for NIR blue SED: absolute magnitude less bright at longer wavelengths, so mag_intrinsic>mag_observed

;;------------------------------------------------------------------
;;-- create grid of observed (unnormalised) magnitudes at different
;;   redshift ranges that we will be probing
;;------------------------------------------------------------------

zgrid = findgen(nbin_kcorr+1)*(maxz_survey-minz_survey)/float(nbin_kcorr)+minz_survey
tmpmag_zgrid = -2.5*alog10(vwsc_conv_filter(wave,spec[*,bestmodel_gal],lam_fil,fil,zgrid))

vmax = fltarr(ngal)

for i=0L,ngal-1 do begin
;   print, i

   if appmag_observed_gal[i] gt maxmag_survey then begin
      print,'VMAX PROBLEM: mag>maxmag_survey, returning -1: galaxy ID ', i
      vmax[i] = -1
      continue
   endif
  if finite(appmag_observed_gal[i]) eq 0 then begin
      print,'VMAX PROBLEM: NAN detected for apparent magnitude, returning -1: galaxy ID ', i
      vmax[i] = -1
      continue
   endif

   if appmag_observed_gal[i] lt minmag_survey then begin
      print,'VMAX PROBLEM: mag<minmag_survey, returning -1: galaxy ID ', i
      vmax[i] = -1
   endif

   if z_gal[i] lt minz_survey then begin
      print,'VMAX PROBLEM: z<minz_survey, returning -1: galaxy ID ', i
      vmax[i] = -1
   endif

   if z_gal[i] gt maxz_survey then begin
      print,'VMAX PROBLEM: z>maxz_survey, returning -1: galaxy ID ', i
      vmax[i] = -1
   endif

;; for the CB: Kcorr to use at grid redshifts 
   tmpmag_zgrid_igal = tmpmag_zgrid[*,i]
   tmpmag_intrinsic_igal = tmpmag_intrinsic[i]
   Kcorr_igal = (tmpmag_intrinsic_igal - tmpmag_zgrid_igal)   

;; K corrections at min z and max z of survey for this galaxy
   Kcorr_max = (tmpmag_intrinsic_igal - tmpmag_zgrid_igal[nbin_kcorr])    
   Kcorr_min = (tmpmag_intrinsic_igal - tmpmag_zgrid_igal[0])    

;;******************************************************************
;;*** upper redshift limit: out to what z can this galaxy be seen? 
;;******************************************************************
;; If the apparent magnitude of this galaxy at zmax is brighter than
;; the survey limit(<maxmag_survey), we don't need to do the
;; calculation. If not, then we have to solve for zmax numerically. 

   if absmag_intrinsic_gal[i]+DM_max-Kcorr_max le maxmag_survey then zmax = maxz_survey else $ ;can be seen at far edge of survey
      if appmag_observed_gal[i] eq maxmag_survey then zmax = z_gal[i] else begin ;galaxy can't be seen at any larger z
      
      zguess = [z_gal[i]+(maxz_survey-z_gal[i])/2.,z_gal[i],maxz_survey]
      appmag = maxmag_survey
      absmag = absmag_intrinsic_gal[i]
      zmax =  fx_root(zguess,'vmax_func',tol=0.01) ;at what z does an absmag galaxy have appmag?

; check for UDSPSB that fx_root is doing something sensible
;      if zmax eq maxz_survey then stop

      ;; fixed by setting higher tolerance
      ;; getting complex solutions? 
      ;; check them and see if ok
      ;; if size(zmax,/type) eq 9 then begin
      ;;    zmax = double(zmax)
      ;;    ind = where(abs(zgrid - zmax) eq min(abs(zgrid - zmax) ))
      ;;    k =  (tmpmag_intrinsic_igal - tmpmag_zgrid_igal[ind])
      ;;    dm = 5.*alog10(lumdist(zmax,/silent)*100000)
      ;;    print,'complex zmax, does this = ', maxmag_survey,'?', i, absmag_intrinsic_gal[i]+DM-k
      ;; endif
   endelse


;;******************************************************************
;;*** lower redshift limit: down to what z can this galaxy be seen? 
;;******************************************************************
;; If the apparent magnitude of this galaxy is fainter than the survey
;; limit (>14.5), we don't need to do the calculation

   if absmag_intrinsic_gal[i]+DM_min-Kcorr_min ge minmag_survey then zmin = minz_survey else $
      if appmag_observed_gal[i] eq minmag_survey then zmin = z_gal[i] else begin

      zguess = [minz_survey+(z_gal[i]-minz_survey)/2.,minz_survey,z_gal[i]]
      appmag = minmag_survey
      absmag = absmag_intrinsic_gal[i]
      zmin = fx_root(zguess,'vmax_func',tol=0.01)

      ;; if size(zmin,/type) eq 9 then begin
      ;;    zmin = double(zmin)
      ;;    ind = where(abs(zgrid - zmin) eq min(abs(zgrid - zmin) ))
      ;;    k =  (tmpmag_intrinsic_igal - tmpmag_zgrid_igal[ind])
      ;;    dm = 5.*alog10(lumdist(zmin,/silent)*100000)
      ;;    print,'complex zmin, does this = ', minmag_survey,'?', i, absmag_intrinsic_gal[i]+DM-k
      ;; endif
   endelse


;; if the galaxy could be seen anywhere in the survey
   if zmin eq minz_survey and zmax eq maxz_survey then begin
      vmax[i] = 1.
      continue
   endif


;; calculate comoving volume corresponding to these redshifts
   vgal = QROMB('comoving',zmin,zmax)

   if vgal lt 0.0 then stop
   
;if vgal/vsurvey lt 0.1 then stop
   
   vmax[i] = vgal/vsurvey
   
endfor

Return, vmax

END
