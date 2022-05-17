;; from Fukugita et al. 1996, eqn 3. (Schneider et al. 1983)

;; see also vw_lambda_eff.pro in same folder

FUNCTION LAMBDA_EFFECTIVE, ll, trans

  c = 2.99792e18                ;Ang/s
  nu = c/ll                     ;Hz

  npix = n_elements(ll)
  dlognu = fltarr(npix) 

  dlognu[0] = alog(nu[1])-alog(nu[0])
  dlognu[npix-1]=alog(nu[npix-1])-alog(nu[npix-2])
  dlognu[1:npix-2]=0.5*(alog(nu[2:npix-1])- alog(nu[0:npix-3]))


  numer=total(trans[0:npix-1]*dlognu*alog(ll),/double)
  denom=total(trans[0:npix-1]*dlognu,/double)
  
  lambda_eff=exp(numer/denom)

  return, lambda_eff
  
END
