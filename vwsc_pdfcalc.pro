;+
; NAME:
;	VWSC_PDFCALC
;
; PURPOSE:
;	This procedure calculates the log-likelihood between a set of
;	models and observations. It builds probability distribution
;	functions, and returns these alongside basic properties of the
;	PDF. 
;
; CATEGORY:
;	Statistics
;
; CALLING SEQUENCE:
;	Write the calling sequence here. Include only positional parameters
;	(i.e., NO KEYWORDS). For procedures, use the form:
;
;	VW_PDFCALC, obs, obs_err, model, params, nbins, pdf_y, pdf_x
;
;
; INPUTS:
;	Obs:	 Observables to be fit.
;       Obs_err: Errors on those observables (ignored if weights are set)
;       Model:   Model values
;       Params:  Structure containing the parameters for which a PDF is required
;       Nbins:   Number of bins to use to calculate the PDF. Array of
;                equal number of elements to params. Can have
;                different number of bins for each parameter to be fit
;
; OPTIONAL INPUTS:
;       Weights: Weights of observables (usually 1/sigma^2 and zero
;                for no data)
;	min:     Minimum value of variable to be fit. Array of
;                equal number of elements to params.
;       max:     Maximum value of variable to be fit. Array of
;                equal number of elements to params.
;
; OUTPUTS:
;	Pdf_y:   structure containing PDF for each variable
;       Pdf_x:   structure containing corresponding variable value for PDF 
;
; OPTIONAL OUTPUTS:
;	stats_pdf:  structure containing mode, 16th, 50th, 84th percentiles of PDF
;       stats_logL: structure containing the index of the
;                   bestfit-model, and best-fit chi-squared value
;
;
;
; EXAMPLE:
;       VWSC_READ_STOCH, dir_model+modeltag, params, spec, lambda,nfile=nfile ;; read in some models
;       ind = (where(params.tlastburst lt 0.5e9 and params.tlastburst gt 0.1e9 and params.fburst[4] gt 0.2))[0:4] ;; select some models to fit
;       flux_mock = flux[*,ind[i]]+randomn(seed,nband)*err[*,ind[i]] ;; add errors to models to be fit
;       VW_PDFCALC, flux_mock, err[*,ind[i]], flux, params_to_fit, nbins,pdf_y, pdf_x, stats_pdf = stats_pdf, stats_logL=stats_logL
;       plot, params[ind].mstr1, lgm_model[1,*],psym=1,xtitle='Input Mass', ytitle = 'Fitted Mass'  ;; compare fitted to input
; 
; MODIFICATION HISTORY:
; 	Written by:	VW 2010 as VW_PDFCALC
;	June, 2014	VW updated to include header and basic checks
;-

PRO VWSC_PDFCALC, obs, obs_err, model, params, nbins, pdf_y, pdf_x, stats_pdf=stats_pdf,stats_logL=stats_logL,weights=weights,min=min,max=max

;;------------------------------------------------------------------
;;-- input parameter checks
;;------------------------------------------------------------------
  nobs = n_elements(obs)
  if nobs eq 1 then begin
     if (size(model,/dim))[0] eq 1 then nmod = (size(model,/dim))[1] else begin
        nmod = (size(model,/dim))[0]
        model = transpose(model)
     endelse
  endif else  nmod = (size(model,/dim))[1]
  
  if nobs gt 1 and nobs ne (size(model,/dim))[0] then stop

  nparam = n_elements(tag_names(params))
  if nmod ne  n_elements(params.(0)) then stop

;;------------------------------------------------------------------
;;-- calculate Likelihood of all models
;;------------------------------------------------------------------
  if not(keyword_set(weights)) then ws = (1/obs_err)^2 else ws = weights
  
  tmp = where(finite(ws) eq 0)
  if tmp[0] ne -1 then ws[tmp]=0.

; slow code 
;  chisq = dblarr(nmod)
;  for j=0,nobs-1 do chisq = chisq+(ws[j]*(reform(model[j,*])-obs[j])^2)

  ;; for testing
;  chisq_mod1 = fltarr(nobs)
;  for j=0,nobs-1 do chisq_mod1[j] = (ws[j]*(reform(model[j,0])-obs[j])^2)
;  splog, 'Chi^2: model 1'
;  splog, chisq_mod1
;  splog, 'total chi^2', total(chisq_mod1)

  
; fast code
  chisq = total(rebin(ws,nobs,nmod)*(model-rebin(obs,nobs,nmod))^2,1,/nan)

  likelihood = exp(-chisq/2.)   ;note that L propto exp(-chisq/2), the constant of proportionality isn't calculated here, 


;;-- output statistics
  stats_logL = {minchisq:min(chisq,/nan),minreducedchisq:min(chisq,/nan)/float(nobs),bestmodelid:(where(chisq eq min(chisq,/nan)))[0]} ;add stats as needed here

;;------------------------------------------------------------------
;;-- create output structure containing the PDF
;;------------------------------------------------------------------
  tagnames = tag_names(params)

  PDF_y = create_struct(tagnames[0],dblarr(nbins[0]))
  for i=1,nparam-1 do PDF_y = create_struct(PDF_y,tagnames[i],dblarr(nbins[i]))

  PDF_x = PDF_y                 ; this is for PDF x-axis array

  stats_pdf = replicate({name:'',per50:0D, mode:0d, per16:0d, per84:0d},nparam)

  for i=0,nparam-1 do begin

     if (size(params.(i)))[0] gt 1 then message, 'One of the parameters is an array - stopping'
     if not(keyword_set(min)) then min_bin = min(params.(i)) else min_bin = min[i]
     if not(keyword_set(max)) then max_bin = max(params.(i)) else max_bin = max[i]

     ;;-- this is the slow bit
     hist = histogram(params.(i),R=R,loc=x,nbin=nbins[i],min=min_bin,max=max_bin)
     xx = x+(x[1]-x[0])/2.

     hist2 = dblarr(nbins[i])
     for j=0,nbins[i]-1 do if R[j] ne R[j+1] then hist2[j] = total(likelihood[R[R[j]:R[j+1]-1]],/nan)

     PDF_y.(i) = hist2         ;the PDF
     PDF_x.(i) = xx          ;the bin positions

     cum_up = total(hist2,/cumulative)
     cum_up = cum_up/max(cum_up) ;0 to 1 in probabiltiy

     stats_pdf[i].name=tagnames[i]
     if max(cum_up) eq 0.0 or (where(finite(cum_up) eq 1))[0] eq -1 then begin
          stats_pdf[i].per50 = -999.
        stats_pdf[i].mode=-999.
        stats_pdf[i].per16 = -999.
        stats_pdf[i].per84=-999.
        continue                ;very bad logL
     endif

;; some useful statistics about the PDF
     tmp = where(hist2 eq max(hist2))
     if n_elements(tmp) eq n_elements(hist2) then stats_pdf[i].mode=-999. $
     else if n_elements(tmp) gt 1 then stats_pdf[i].mode=mean(xx[tmp]) $
     else stats_pdf[i].mode = xx[tmp]

     stats_pdf[i].per50 = xx[(where(cum_up gt 0.5))[0]]
     stats_pdf[i].per16 = xx[(where(cum_up gt 0.16))[0]]
     stats_pdf[i].per84 = xx[(where(cum_up gt 0.84))[0]]
     
  endfor

END
