FUNCTION combinesed, evecfile, zbin, flux_super_data_norm,  flux_super_data_err_norm, zphot

  nzbin = n_elements(zbin)-1
  npix = (size(flux_super_data_norm,/dim))[0]

  restore, evecfile ; need evecs for projection of avspec
  namplitudes = (size(pcs,/dim))[0]

;;-- set up average spectrum structure
  avspec = replicate({sed:fltarr(npix),err:fltarr(npix),pcs:fltarr(namplitudes), pcerr:fltarr(namplitudes)},nzbin) 

  for j=0,nzbin-1 do begin
     ind_z = where(zphot gt zbin[j] and zphot lt zbin[j+1],nn)

;;-- calculate 16th, 50th and 84th percentiles and overplot
     for i=0,npix-1 do begin
        tmp = where(flux_super_data_norm[i,ind_z] gt 0,n)
        if n eq 1 then begin ;one object in this lambda bin
           avspec[j].sed[i] = flux_super_data_norm[i,ind_z[tmp]]
           avspec[j].err[i] = flux_super_data_err_norm[i,ind_z[tmp]]
        endif
        if n gt 1 then begin
           avspec[j].sed[i] = mean(flux_super_data_norm[i,ind_z[tmp]])
           avspec[j].err[i] = sqrt(total(flux_super_data_err_norm[i,ind_z[tmp]]^2))/float(n)
        endif
     endfor
  
;;-- calculate PCs of this average SED, using percentiles as 1sigma errors  
     avspec[j].pcs = vwpca_normgappy(avspec[j].sed, avspec[j].err, evecs[*,0:namplitudes-1], meanarr, cov=cov)
     avspec[j].pcerr = sqrt(diag_matrix(cov[*,*]))

     print, avspec[j].pcs
     print, avspec[j].pcerr
     ploterror, wave_rest_super, avspec[j].sed, avspec[j].err,psym=1,title=string(zbin[j],form='(F0.1)')+'<z<'+string(zbin[j+1],form='(F0.1)')

  endfor
  

return, avspec
END

;;------------------------------------------------------------------
;;------------------------------------------------------------------

PRO VWSC_FIGS_PDF, dir_figs, evecfile, superfluxfile, norm_data, pcs_data, pcerr_data, lines, pcs_models, params_models, zphot, minz,maxz,nzbin

  restore, superfluxfile ;contains flux_super_data (and errors, but not used here)

  zbin = fltarr(nzbin+1)
  for i=0,nzbin do zbin[i] = i*(maxz-minz)/float(nzbin)+minz ;end of last bin

;;-- useful sizes
  nmod = n_elements(params_models)
  ngal = (size(flux_super_data,/dim))[1]
  npix = (size(flux_super_data,/dim))[0]

;;-- Classes
  VWSC_CLASSIFY, pcs_data, lines,ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty, ind_lomet,ind_junk, class, stack_name, stack_color

;;-- normalise flux arrays, using PCA fitted flux normalisation
  flux_super_data_norm = flux_super_data/transpose(rebin(norm_data,ngal,npix))
  flux_super_data_err_norm = flux_super_data_err/transpose(rebin(norm_data,ngal,npix))

;;-- Calculate average PSB SED and best-fit PCs
  ps1, dir_figs+'SED_combine_psb.ps'
  avspec_psb = combinesed(evecfile,zbin,flux_super_data_norm[*,ind_psb],  flux_super_data_err_norm[*,ind_psb], zphot[ind_psb])
  ps2

;;-- Calculate and plot PDF 
  params_fit = replicate({tform:0d, gamma:0d, zmet:0d, tauv0:0d, tlastburst:0d,fburst_gyr:0d,age_wm:0d, age_wr:0d},nmod)
  copy_struct, params_models, params_fit
  params_fit.fburst_gyr = params_models.fburst[3]
  params_fit.tlastburst /= 1d9
  params_fit.age_wr /= 1e9
  params_fit.age_wm /= 1e9
  params_fit.tform /= 1e9

  ps1, dir_figs+'PDF_params_psb.ps'
  ind = [0,1,2,3,4,5,6,7] ;[1,2,3,4,5,9,15]
  name = tag_names(params_fit)

  nparam = n_elements(name)
  nbins_pdf=fltarr(nparam)+10
  maxbin = fltarr(nparam)
  for i=0,nparam-1 do maxbin[i] = max(params_fit.(i))
  maxbin[where(name eq 'TLASTBURST')] = 2.
  maxbin[where(name eq 'AGE_WR')] = 8.
  maxbin[where(name eq 'TAUV0')] = 5.

  device,/portrait,xoffset=0,yoffset=0,ysize=20,xsize=25
  !p.multi=[0,n_elements(ind),nzbin]
  !x.margin= [4,3]
  !y.margin= [2,2]
  !y.omargin=[2,0]
  !x.omargin=[7,0]
 
  for j=0,nzbin-1 do begin

     vw_pdfcalc,avspec_psb[j].pcs,avspec_psb[j].pcerr,pcs_models,params_fit,nbins_pdf,pdf_psb,xx_psb,max=maxbin ;,stats_pdf=stats_pdf,stats_logL=stats_logL ;,min=min,max=max)
     
     for k=0,n_elements(ind)-1 do plot, xx_psb.(ind[k]), pdf_psb.(ind[k])/total(pdf_psb.(ind[k])),xtitle=name[ind[k]],psym=10,/xs
     xyouts,0.025, (nzbin-j)/float(nzbin)*0.95, string(zbin[j],form='(F0.1)')+'<z<'+string(zbin[j+1],form='(F0.1)'),/norm,align=0,orien=90

  endfor

  ps2

END

;;-- this didn't appear to work, but now see I was fitting
;;   everything, not just PSBs!

;;   params_fit = replicate({tform:0d, gamma:0d, zmet:0d, tauv0:0d, tlastburst:0d,fburst_gyr:0d,age_wm:0d, age_wr:0d},nmod)
;;   copy_struct, params, params_fit
;;   params_fit.fburst_gyr = params.fburst[3]
;;   params_fit.tlastburst /= 1d9
;;   params_fit.age_wr /= 1e9
;;   params_fit.age_wm /= 1e9
;;   params_fit.tform /= 1e9

;;   ps1, dir_figs+'PDF_params.ps'
;;   ind = [0,1,2,3,4,5,6,7];[1,2,3,4,5,9,15]
;;   name = tag_names(params_fit)

;;   nparam = n_elements(name)
;;   nbins_pdf=fltarr(nparam)+20

;;   device,/portrait,xoffset=0,yoffset=0,ysize=20,xsize=25
;;   !p.multi=[0,n_elements(ind),nzbin]
;;   !x.margin= [4,3]
;;   !y.margin= [2,2]
;;   !y.omargin=[2,0]
;;   !x.omargin=[7,0]

;; ;  multiplot,[7,7],mxtitle=textoidl('Rest wavelength,\mum'),mytitle=textoidl('Normalised Flux'),ygap=0.01,xgap=0.01,mytitoffset=-1,/doyaxis
;;   for j=0,nzbin-1 do begin
;;      ind_z = where(zphot gt zbin[j] and zphot lt zbin[j+1],nn)
;;      print, 'combining PDFs for N galaxies in zbin:', nn,zbin[j],zbin[j+1]
;;      for i=0,nn-1 do begin
;;         ;;  obs, obs_err, model, params, nbins, pdf_y, pdf_x, 
;;         vw_pdfcalc,pcs_data[*,ind_z[i]],pcerr_data[*,ind_z[i]],pcs,params_fit,nbins_pdf,pdf,xx;,stats_pdf=stats_pdf,stats_logL=stats_logL ;,min=min,max=max)
;;         if i eq 0 then begin
;;            xx_psb = xx          ;should be same for all galaxies
;;            pdf_psb = pdf        ;initialise PDF array
;;            for k=0,n_elements(ind)-1 do pdf_psb.(ind[k]) = 1.
;;         endif

;;-- this is not correct - should multiply probabilities
;;   http://www.mpia.de/Gaia/publications/probcomb_TN.pdf eqn. 8 USEFUL/CombiningProbabilities.pdf
;;         for k=0,n_elements(ind)-1 do pdf_psb.(ind[k]) = pdf_psb.(ind[k])+pdf.(ind[k])/total(pdf.(ind[k]))   ;multiply all galaxies together
;;         ;plot, xx_psb.(ind[0]), pdf_psb.(ind[0]) ,psym=10
;;      endfor


;;      for k=0,n_elements(ind)-1 do plot, xx_psb.(ind[k]), pdf_psb.(ind[k])/total(pdf_psb.(ind[k])),xtitle=name[ind[k]],psym=10,/xs
;;      xyouts,0.025, (nzbin-j)/float(nzbin)*0.9, string(zbin[j],form='(F0.1)')+'<z<'+string(zbin[j+1],form='(F0.1)'),/norm,align=0,orien=90

  


;;   endfor

;;   ps2

