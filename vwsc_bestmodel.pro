FUNCTION VWSC_BESTMODEL,  obs, obs_err, model, minchisq=minchisq

;;-- check only got one object
  if (size(obs))[0] gt 1 then begin
     print, 'VWSC_BESTMODEL: please input one object at a time'
     return, -1
  endif

  if obs_err[0] eq 0.0 then begin
      print, 'VWSC_BESTMODEL: zero error, returning -1'
      minchisq=0.0
      return,-1
  endif

;;-- numbers and checks
  nobs = n_elements(obs)
  if nobs eq 1 then begin
     if (size(model,/dim))[0] eq 1 then nmod = (size(model,/dim))[1] else begin
        nmod = (size(model,/dim))[0]
        model = transpose(model)
     endelse
  endif else  nmod = (size(model,/dim))[1]
  
  if nobs gt 1 and nobs ne (size(model,/dim))[0] then stop

;;-- calculate chi^2 of all models relative to this observation
  chisq = dblarr(nmod)
  for j=0,nobs-1 do if obs_err[j] ne 0 then chisq = chisq+((reform(model[j,*])-obs[j])/obs_err[j])^2

  chisq = chisq/float(nobs)     ;reduced chisq

  minchisq = min(chisq,/nan)
  
  if minchisq eq 0.0 then bestmodelid = -1

  bestmodelid = where(chisq eq min(chisq,/nan))

  RETURN, bestmodelid

END
