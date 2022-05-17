
FUNCTION VWSC_PROJECTDATA, gridfile, evecfile, filterset, dir_figs,ID, flux, flux_err, z, namplitudes, minz, maxz, dz, $
                           chisq=chisq,ngood=ngood,pcerr=pcerr,norm=norm,outfile=outfile,superfluxfile=superfluxfile, $
                           maxplot=maxplot

  restore, gridfile            
  restore, evecfile
  npix = n_elements(wave_rest_super)


;;-- filters
  readcol, filterset, filterfiles,ll_eff,filternames,form='(A, F,A)',/silent ;list of filter files
  nband = n_elements(ll_eff)


  ngal = (size(flux,/dim))[1]
  splog, '# galaxies: ',ngal

;  if ngal ge 100 then  maxplot = 100 else maxplot = ngal               ;maximum number of objects to plot

  
  if  (size(flux,/dim))[0] ne nband then stop
  
;;-- place f_nu_obs into super-sampled array, convert into f_lambda_rest
  flux_super_data = (flux_super_data_err = fltarr(npix,ngal))
  for i=0l,ngal-1 do flux_super_data[*,i] = vwsc_fillflux(flux[*,i],z[i],minz,maxz,dz,ll_eff,ind_wave,/jansky,/observed)
  for i=0l,ngal-1 do flux_super_data_err[*,i] = vwsc_fillflux(flux_err[*,i],z[i],minz,maxz,dz,ll_eff,ind_wave,/jansky,/observed)
  
;;-- f_lambda_rest
  
  if keyword_set(superfluxfile) then save, flux_super_data, flux_super_data_err,wave_rest_super,file=superfluxfile

  project_data = uintarr(npix,ngal)
  project_data[where(flux_super_data_err ne 0)] = 1

;;-- evecs are in f_lambda_rest -> norm_data is in f_lambda_rest 
;;-- normalise close to 1 for projection to prevent projection errors

  pcs_data = vwpca_normgappy(flux_super_data/1d5,flux_super_data_err/1d5, evecs[*,0:namplitudes-1], meanarr,norm=norm,cov=cov)
  norm = norm*1d5
  pcerr = fltarr(namplitudes,ngal)
  for i=0L,ngal-1 do pcerr[*,i] = sqrt(diag_matrix(cov[*,*,i]))


;;-- calculate chisq of projection
  chisq = fltarr(ngal)
  ngood = fltarr(ngal)
  for i=0L,ngal-1 do begin
     
     ind_good = where(project_data[*,i] eq 1,nn)
     ngood[i] = nn              ;  # of bands
     
     recon = (vwpca_reconstruct(pcs_data[*,i],evecs[*,0:namplitudes-1])+meanarr)*norm[i]  
     chisq[i] = total((recon[ind_good]-flux_super_data[ind_good,i])^2/(flux_super_data_err[ind_good,i])^2) / float(nn)

  endfor

;;-- example spectra
if maxplot gt 0 then begin
  ps1c,dir_figs+'vwsc_projectdata.ps'
  !p.multi=[0,1,5]

  for i=0,maxplot-1 do begin
     recon = (vwpca_reconstruct(pcs_data[*,i],evecs[*,0:namplitudes-1])+meanarr)*norm[i]
     ind = where(flux_super_data_err[*,i] gt 0. and flux_super_data_err[*,i] lt 10000000.)
        
     ploterror, wave_rest_super[ind], flux_super_data[ind,i],flux_super_data_err[ind,i],xrange=minmax(wave_rest_super),/xs,/xlog,psym=2,title='ID='+string(ID[i],form='(I0.0)')+' pcs='+string(pcs_data[*,i],form='('+string(namplitudes,form='(I0)')+'(F0.1,1X))')+' pcerr='+string(pcerr[*,i],form='('+string(namplitudes,form='(I0)')+'(F0.1,1X))')+' z='+string(z[i],form='(F0.2)')+' chisq='+string(chisq[i],form='(F0.2)')
     
     oplot, wave_rest_super,recon,psym=3
     oplot, wave_rest_super[ind], recon[ind],psym=4
     if n_elements(filternames_super) ne 0 then xyouts, wave_rest_super[ind],replicate(median(recon[ind]),nband),filternames_super[ind]
     
  endfor 

  !p.multi=[0,1,2]
  plot_histogram, chisq,nbins=100,max=20 , title='chisq' 
  ps2

endif     
 
if n_elements(outfile) ne 0 then begin
   openw,1,outfile,width=200
;   printf,1,'# ID mass ra dec pc1 pc2 pc3 pc1err pc2err pc3err chisq z nband'
;   for i=0,ngal-1 do printf,1,ID[i], mass[i], ra[i],dec[i],pcs_data[0,i],pcs_data[1,i],pcs_data[2,i],pcerr[0,i],pcerr[1,i],pcerr[2,i],chisq[i],z[i],ngood[i],form='(I0,1X,F0.2,1X,2(F0.6,1X),9(F0.2,1X))'
   printf,1,'# ID pc1 pc2 pc3 pc1err pc2err pc3err chisq z nband'
   for i=0,ngal-1 do printf,1,ID[i],pcs_data[0,i],pcs_data[1,i],pcs_data[2,i],pcerr[0,i],pcerr[1,i],pcerr[2,i],chisq[i],z[i],ngood[i],form='(I0,1X,9(F0.2,1X))'

   close,1
endif

return, pcs_data
END
