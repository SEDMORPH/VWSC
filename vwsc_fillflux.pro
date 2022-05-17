;;*** AIM: take observed flux arr and redshift, and place into
;;         correct indices of supersampled flux arr for PCA
;; 
;;*** INPUT: observed flux of single object in f_lambda [or f_nu if /jansky keyword is set]
;;           ll_obs : observed wavelength of observations in AA
;;           z : redshift
;;
;;*** OUTPUT: rest-frame fluxarr (f_lambda_rest)
;;
FUNCTION VWSC_FILLFLUX, flux,z,minz,maxz,dz,ll_obs,ind_wave,jansky=jansky,observed=observed

  @vwsc_constants.inc

  nredshift = round((maxz-minz)/dz)+1
  zbin = findgen(nredshift)*dz + minz 

  nz = n_elements(zbin)
  nband = n_elements(ll_obs)

;;-- f_nu_obs to f_lambda_obs if required
  if keyword_Set(jansky) then  ff = c_in_aa*flux/ll_obs^2 else ff = flux
  
;;-- rest-frame
  if keyword_set(observed) then ff = ff*(1+z)
  
;;-- which z bin? 
  tmp = abs(z-zbin)
  ind_zz = (where(tmp eq min(tmp)))[0]

;;-- which bins to put flux into 
  ind = findgen(nband)*nz+ind_zz

  fluxarr = fltarr(nz*nband)
  fluxarr[ind] = ff
  
  fluxarr  = fluxarr[ind_wave]

  return, fluxarr

END
