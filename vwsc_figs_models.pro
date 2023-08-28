
PRO PLOT_PARAMS, psfile, x,y,params,ABmag,xtitle,ytitle,range,bins,_extra=extra,loadct=loadct,ctreverse=ctreverse

ps1,psfile
   device,/portrait,xoffset=0,yoffset=0,ysize=19.5,xsize=20
ncolors=100
loadct,10
plotsym,0,0.1,/fill
!x.margin=[8,3]
multiplot,[2,3]

mag_sun = 4.55

mass = alog10(params.mstr1*10.^(0.4*(ABmag-mag_sun)))-0.4*(ABmag-mag_sun)

ind = where(params.age_wr ge 7.2 and params.zmet le 2 and params.zmet ge 0 and params.tauv0 ge 0 and params.tauv0 le 6 and params.fburst[3] ge 0 $
and params.fburst[3] le 1 and x le 165)
params_clean = params[ind]
mass_clean = mass[ind]
x_clean = x[ind]
y_clean = y[ind]

nmod = n_elements(params_clean)

vwplot_2dhist,x_clean,y_clean,range,bins, weight1=alog10(params_clean.age_wr),vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out, ytitle=ytitle,ctreverse=ctreverse
;color = bytscl(params.age_wr,top=ncolors-1)
;plot, x,y,/nodata,ytitle=ytitle,_extra=extra;,position = [0.15,0.55,0.45,0.95]
;plots,x,y,color=color,psym=8
;cgcolorbar,position=[0.46,0.55,0.48,0.90],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right

;cgcolorbar,position=[0.46,0.55,0.48,0.85],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
cgcolorbar,position=[0.46,0.7,0.48,0.9],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
xyouts, 0.47,0.92,textoidl('log(Age_{wr})'),/norm,align=0.5,charsize=0.8

multiplot
vwplot_2dhist,x_clean,y_clean,range,bins, weight1=params_clean.zmet,vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out,ctreverse=ctreverse
;color = bytscl(params.zmet,top=ncolors-1)
;plot, x,y,/nodata,_extra=extra;,position = [0.55,0.55,0.85,0.95]
;plots,x,y,color=color,psym=8
cgcolorbar,position= [0.87,0.7,0.89,0.9],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
xyouts, 0.9,0.92,textoidl('Z'),/norm,align=0.5,charsize=0.8

multiplot
vwplot_2dhist,x_clean,y_clean,range,bins, weight1=params_clean.tauv0,vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out, ytitle=ytitle,ctreverse=ctreverse
;color = bytscl(params.tauv0,top=ncolors-1)
;plot, x,y,/nodata,xtitle=xtitle,ytitle=ytitle,_extra=extra;,position = [0.15,0.15,0.45,0.5]
;plots,x,y,color=color,psym=8

;cgcolorbar,position= [0.46,0.15,0.48,0.46],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
cgcolorbar,position= [0.46,0.43,0.48,0.63],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
xyouts, 0.47,0.65,textoidl('\tau_v'),/norm,align=0.5,charsize=0.8

multiplot,/doxaxis
vwplot_2dhist,x_clean,y_clean,range,bins, weight1=params_clean.fburst[3],vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out,xtitle=xtitle
;;params.tlastburst/1e9
;color = bytscl(params.tlastburst,min=0,max=1e9,top=ncolors-1)
;plot, x,y,/nodata,xtitle=xtitle,_extra=extra;,position = [0.55,0.15,0.85,0.5]
;plots,x,y,color=color,psym=8
;; time: 1.e6, 1.e7, 1.e8, 1.e9 and 2.e9
;ind = where(params.fburst[3] gt 0.2 and params.tlastburst lt 1e9)
;plots,x[ind],y[ind],color=color[ind],psym=1

cgcolorbar,position= [0.87,0.43,0.89,0.63],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.82,divisions=5,/reverse,/right
xyouts, 0.9,0.65,textoidl('F_{burst}'),/norm,align=0.5,charsize=0.8

;multiplot
;vwplot_2dhist,x,y,range,bins, weight1=mass_clean,vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out, ytitle=ytitle,xtitle=xtitle,ctreverse=ctreverse
;color = bytscl(params.tauv0,top=ncolors-1)
;plot, x,y,/nodata,xtitle=xtitle,ytitle=ytitle,_extra=extra;,position = [0.15,0.15,0.45,0.5]
;plots,x,y,color=color,psym=8

;cgcolorbar,position= [0.46,0.15,0.48,0.35],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
;xyouts, 0.47,0.37,textoidl('log(M*/M_{solar})'),/norm,align=0.5,charsize=0.8

;The next lines commented out because the rand files don't contain sfrav, so the code crashes. Plot not used so commented out.
;multiplot
;vwplot_2dhist,x_clean,y_clean,range,bins, weight1=alog10(params_clean.sfrav[2]/params_clean.mstr1),vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out,xtitle=xtitle,ytitle=ytitle,ctreverse=ctreverse
;color = bytscl(params.tauv0,top=ncolors-1)
;plot, x,y,/nodata,xtitle=xtitle,ytitle=ytitle,_extra=extra;,position = [0.15,0.15,0.45,0.5]
;plots,x,y,color=color,psym=8

cgcolorbar,position= [0.46,0.15,0.48,0.35],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
xyouts, 0.47,0.37,textoidl('log(sSFR)'),/norm,align=0.5,charsize=0.8




multiplot,/reset


ps2

END


;;****************************************************************************************
;;****************************************************************************************
PRO VWSC_FIGS_MODELS, dir_figs, dir_notes, gridfile, mockfile, evecfile, filterset, lambda, spec, params, namplitudes,wave_norm,plotparams=plotparams
  
  restore, gridfile
  restore, mockfile
  restore, evecfile
  
  wavemin = min(wave_rest_super)
  wavemax = max(wave_rest_super)

  splog, file=dir_notes+'vwsc_figs_models.splog'

;;-- filters
  readcol, filterset, filterfiles,ll_eff,filternames,form='(A, F,A)',/silent ;list of filter files
  nband = n_elements(ll_eff)
 
;;------------------------------------------------------------------
;;-- Figures
;;------------------------------------------------------------------

;;-- e.g. spectra

  ps1,dir_figs+'mockspec.ps'
  plotsym,0,0.15,/fill
  !x.margin=[6,3]
  multiplot,[1,5]

  plotnames = strmid(filternames,0,1)
  ind = where(plotnames eq '3') & if ind[0] ne -1 then plotnames[ind] = '3.6'
  ind = where(plotnames eq '4') & if ind[0] ne -1 then plotnames[ind] = '4.5'

  ind_plot = [19,13,0,1,4]        ;random selection
  
  splog, 'redshifts of eg spectra: ', zmod_phot[ind_plot]

  for i=0,4 do begin
     if i ne 0 then multiplot
     if i eq 4 then xtitle=textoidl('Rest-frame wavelength, /\mum') else xtitle = ''
     if i eq 2 then ytitle=textoidl('Normalised F_\lambda') else ytitle = ''
     
;        yr = [min(flux_super[*,i])*0.9,max(flux_super[*,i])*1.1] 
     yr = [0,max(flux_super[*,ind_plot[i]])*1.1]
     if yr[1] gt 2 and yr[1] lt 2.1 then yr[1]=2.3
     
     if max(yr) lt 4 then ytickint=1 else ytickint = 2
     plot, wave_rest_super/1e4,flux_super[*,ind_plot[i]], psym=8,xtitle=xtitle,ytitle=ytitle,$
           xr=[wavemin-500,wavemax+100]/1e4,/xs,ytickinterval=ytickint,yr=yr,/ys
     
     ;; add diamonds of real filters
;     ind = where(ll_eff/(1+zmod[ind_plot[i]]) gt wavemin and ll_eff/(1+zmod[ind_plot[i]]) lt wavemax and flux_uds[*,ind_plot[i]] ne 0,nn)
;     oplot, ll_eff[ind]/(1+zmod[ind_plot[i]])/1e4,flux_uds[ind,ind_plot[i]],psym=4
     
     ;; add filternames, except 3.6
;     ind = where(ll_eff/(1+zmod[ind_plot[i]]) gt wavemin and ll_eff/(1+zmod[ind_plot[i]]) lt wavemax and plotnames ne '3.6' and flux_uds[*,ind_plot[i]] ne 0,nn)
;     xyouts, ll_eff[ind]/(1+zmod[ind_plot[i]])/1e4,fltarr(nn)+(min(yr)+(max(yr)-min(yr))*0.05),plotnames[ind],align=0.5,color=cgcolor('dark gray'),charsize=0.7
     
     ;; add 3.6 band filtername, slightly offset
;     ind = where(ll_eff/(1+zmod[ind_plot[i]]) gt wavemin and ll_eff/(1+zmod[ind_plot[i]]) lt wavemax and plotnames eq '3.6' and flux_uds[*,ind_plot[i]] ne 0)
;     if ind[0] ne -1 then xyouts, ll_eff[ind]/(1+zmod[ind_plot[i]])/1e4-0.03,(min(yr)+(max(yr)-min(yr))*0.05),plotnames[ind],align=0.5,color=cgcolor('dark gray'),charsize=0.7
     
;     xyouts, (wavemax-500)/1e4, min(yr)+(max(yr)-min(yr))*0.8,'z = '+string(zmod[ind_plot[i]],form='(f0.1)'),align=1
     


  endfor
 
  multiplot,/reset
  ps2


;;-- E-vectors
  ps1,dir_figs+'evecs.ps'
  xx = 0.9
  xr = [wavemin-500,wavemax+100]/1e4
  !y.charsize=0.9
  ;!x.margin=[6,3]
  multiplot, [1,namplitudes+1]
  plot, wave_rest_super/1e4,meanarr,psym=8,xr=xr,/xs;,yr=[0,2.5],/ys
  xyouts, xx,0.9,'Mean',/norm,align=1
  multiplot
  plot, wave_rest_super/1e4,evecs[*,0],psym=8,xr=xr,/xs,/ys,ytickinterval=0.02;,yr=[-0.01,0.07]
  xyouts, xx,0.68,'E-vector 1',/norm,align=1
  multiplot
  plot, wave_rest_super/1e4,evecs[*,1],psym=8,xr=xr,/xs,/ys,ytickinterval=0.04;,yr=[-0.09,0.05]
  xyouts, xx,0.47,'E-vector 2',/norm,align=1

  multiplot
  plot, wave_rest_super/1e4,evecs[*,2],psym=8,xr=xr,/xs,xtitle=textoidl('Rest-frame wavelength, /\mum'),/ys,ytickinterval=0.04;,yr=[-0.09,0.06]
  xyouts, xx,0.26,'E-vector 3',/norm,align=1
  
  xyouts, 0.05,0.5,textoidl('Normalised F_\lambda'),orien=90,/norm,align=0.5

  multiplot,/reset
  ps2
  

;;-- Delta PCs
  ps1,dir_figs+'deltapcs.ps'
  min = -1
  max = 1
  yr = [0,0.15]
  
  multiplot,[namplitudes,1]
  plot_histogram, pcs_mock[0,*]-pcs[0,*],min=min,max=max,nbins=50,linestyle=1,xtitle=textoidl('\Delta(SC1)'),ytitle='Fraction of mock galaxies',yr=yr,/ys,/norm
  plot_histogram, pcs_mock_zphot[0,*]-pcs[0,*],min=min,max=max,nbins=50,/over,linestyle=2,/norm
  plots, [0.15,0.18],[0.9,0.9],/norm & xyouts, 0.2,0.9,'(a)',align=0,/norm
  plots, [0.15,0.18],[0.85,0.85],/norm,linestyle=1 & xyouts, 0.2,0.85,'(b)',align=0,/norm
  plots, [0.15,0.18],[0.8,0.8],/norm,linestyle=2 & xyouts, 0.2,0.8,'(c)',align=0,/norm
  
  multiplot
  plot_histogram, pcs_mock[1,*]-pcs[1,*],min=min,max=max,nbins=50,linestyle=1,/norm,xtitle=textoidl('\Delta(SC2)'),yr=yr,/ys
  plot_histogram, pcs_mock_zphot[1,*]-pcs[1,*],min=min,max=max,nbins=50,/over,linestyle=2,/norm

  multiplot
  plot_histogram, pcs_mock[2,*]-pcs[2,*],min=min,max=max,nbins=50,linestyle=1,/norm,xtitle=textoidl('\Delta(SC3)'),yr=yr,/ys
  plot_histogram, pcs_mock_zphot[2,*]-pcs[2,*],min=min,max=max,nbins=50,/over,linestyle=2,/norm
  
  multiplot,/reset
  ps2

  splog, 'SC1 percentiles:', vw_percentile(pcs_mock_zphot[0,*]-pcs[0,*])
  splog, 'SC2 percentiles:', vw_percentile(pcs_mock_zphot[1,*]-pcs[1,*])
  splog, 'SC3 percentiles:', vw_percentile(pcs_mock_zphot[2,*]-pcs[2,*])
  ind_loz = where(zmod_phot gt 0.9 and zmod_phot lt 1.2)
  splog, 'SC1 loz percentiles:', vw_percentile(pcs_mock_zphot[0,ind_loz]-pcs[0,ind_loz])
  splog, 'SC2 loz percentiles:', vw_percentile(pcs_mock_zphot[1,ind_loz]-pcs[1,ind_loz])
  splog, 'SC3 loz percentiles:', vw_percentile(pcs_mock_zphot[2,ind_loz]-pcs[2,ind_loz])
  splog, 'SC1 loz relative error percentiles:', vw_percentile((pcs_mock_zphot[0,ind_loz]-pcs[0,ind_loz])/pcs[0,ind_loz])
  splog, 'SC2 loz relative error percentiles:', vw_percentile((pcs_mock_zphot[1,ind_loz]-pcs[1,ind_loz])/pcs[1,ind_loz])
  splog, 'SC3 loz relative error  percentiles:', vw_percentile((pcs_mock_zphot[2,ind_loz]-pcs[2,ind_loz])/pcs[2,ind_loz])

;;-- redshift dependence of information
  ps1,dir_figs+'deltapcs_z.ps'
  multiplot,[1,3]
  dpc = pcs_mock_zphot[0,*]-pcs[0,*]
  plot, zmod,dpc,psym=3,ytitle=textoidl('\Delta(SC1)'),/xs,yr=[-1.54,1.54],/ys
  hist = histogram(zmod,binsize=0.05,reverse_ind=R,loc=x) & nbin = n_elements(hist)
  per = fltarr(nbin,3)
  for i=0,nbin-1 do IF R[i] NE R[i+1] THEN per[i,*] = vw_percentile(dpc[R[R[I] : R[i+1]-1]])
  oplot, x+(x[1]-x[0])/2.,per[*,0],color=cgcolor('cyan'),linestyle=2,thick=6
  oplot, x+(x[1]-x[0])/2.,per[*,1],color=cgcolor('cyan'),thick=6
  oplot, x+(x[1]-x[0])/2.,per[*,2],color=cgcolor('cyan'),linestyle=2,thick=6
  multiplot
  
  dpc = pcs_mock_zphot[1,*]-pcs[1,*]
  plot, zmod,dpc,psym=3,ytitle=textoidl('\Delta(SC2)'),/xs,yr=[-1.54,1.54],/ys
  per = fltarr(nbin,3)
  for i=0,nbin-1 do IF R[i] NE R[i+1] THEN per[i,*] = vw_percentile(dpc[R[R[I] : R[i+1]-1]])
  oplot, x+(x[1]-x[0])/2.,per[*,0],color=cgcolor('cyan'),linestyle=2,thick=6
  oplot, x+(x[1]-x[0])/2.,per[*,1],color=cgcolor('cyan'),thick=6
  oplot, x+(x[1]-x[0])/2.,per[*,2],color=cgcolor('cyan'),linestyle=2,thick=6
  multiplot

  dpc = pcs_mock_zphot[2,*]-pcs[2,*] 
  plot, zmod,dpc,psym=3,xtitle='z',ytitle=textoidl('\Delta(SC3)'),/xs,yr=[-1.54,1.54],/ys
  per = fltarr(nbin,3)
  for i=0,nbin-1 do IF R[i] NE R[i+1] THEN per[i,*] = vw_percentile(dpc[R[R[I] : R[i+1]-1]])
  oplot, x+(x[1]-x[0])/2.,per[*,0],color=cgcolor('cyan'),linestyle=2,thick=6
  oplot, x+(x[1]-x[0])/2.,per[*,1],color=cgcolor('cyan'),thick=6
  oplot, x+(x[1]-x[0])/2.,per[*,2],color=cgcolor('cyan'),linestyle=2,thick=6
  multiplot,/reset
  ps2

;;-- reconstructions + residuals


  ps1,dir_figs+'mockrecon.ps'
  plotsym,0,0.15,/fill
;  !x.margin=[6,3]
  multiplot,[1,5]
  for i=0,4 do begin
     if i ne 0 then multiplot
     if i eq 4 then xtitle=textoidl('Rest-frame wavelength, /\mum') else xtitle = ''

     yr = [0,max(flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]])*1.1]
     if yr[1] gt 2 and yr[1] lt 2.1 then yr[1]=2.3
     if max(yr) lt 4 then ytickint=1 else ytickint = 2

     plot, wave_rest_super/1e4,flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]], psym=8,xtitle=xtitle,ytitle=ytitle,$
           xr=[wavemin-500,wavemax+100]/1e4,/xs,yrange=yr,/ys,ytickinterval=ytickint
     
     recon = (vwpca_reconstruct(pcs_mock_zphot[*,ind_plot[i]],evecs[*,0:namplitudes-1])+meanarr)
     oplot, wave_rest_super/1e4,recon,psym=8,color=cgcolor('green')

    ;; add diamonds of real filters
     ind = where(ll_eff/(1+zmod[ind_plot[i]]) gt wavemin and ll_eff/(1+zmod[ind_plot[i]]) lt wavemax and flux_mock[*,ind_plot[i]] ne 0,nn)
     oplot, ll_eff[ind]/(1+zmod[ind_plot[i]])/1e4,flux_mock[ind,ind_plot[i]]/norm_super[ind_plot[i]],psym=4
 
     xyouts, (wavemax-500)/1e4, min(yr)+(max(yr)-min(yr))*0.8,'z = '+string(zmod[ind_plot[i]],form='(f0.1)'),align=1

 
  endfor
  xyouts, 0.05,0.5,textoidl('Normalised F_\lambda'),/norm,align=0.5,orien=90

  ps2


  ps1,dir_figs+'mockresid.ps'
  plotsym,0,0.15,/fill
;  !x.margin=[6,3]
  multiplot,[1,5]
  for i=0,4 do begin
     if i ne 0 then multiplot
     if i eq 4 then xtitle=textoidl('Rest-frame wavelength, /\mum') else xtitle = ''

     yr = [-0.05,0.045]

     recon = (vwpca_reconstruct(pcs[*,ind_plot[i]],evecs[*,0:namplitudes-1])+meanarr)
     plot, wave_rest_super/1e4,(flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]]-recon)/(flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]]), psym=8,xtitle=xtitle,ytitle=ytitle,ytickint=0.05, $
           xr=[wavemin-500,wavemax+100]/1e4,/xs,yr=yr,/ys;,yrange=[min(flux_super[*,ind_plot[i]])*0.9,max(flux_super[*,ind_plot[i]])*1.3],/ys,ytickinterval=1
     
     recon = (vwpca_reconstruct(pcs_mock[*,ind_plot[i]],evecs[*,0:namplitudes-1])+meanarr)
     oplot, wave_rest_super/1e4,(flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]]-recon)/(flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]]),psym=8,color=cgcolor('blue')

     recon = (vwpca_reconstruct(pcs_mock_zphot[*,ind_plot[i]],evecs[*,0:namplitudes-1])+meanarr)
     oplot, wave_rest_super/1e4,(flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]]-recon)/(flux_super[*,ind_plot[i]]/norm_super[ind_plot[i]]),psym=8,color=cgcolor('red')

;     ind = where(flux_super_uds[*,ind_plot[i]] ne 0)
;     print, recon[ind]-flux_super_uds[ind,ind_plot[i]]

     ;; add diamonds of real filters
     ind = where(flux_super_mock[*,ind_plot[i]] ne 0,nn)
     oplot, wave_rest_super[ind]/1e4,(flux_super_mock[ind,ind_plot[i]]/norm_mock[ind_plot[i]]-recon[ind])/(flux_super_mock[ind,ind_plot[i]]/norm_mock[ind_plot[i]]),psym=4,color=cgcolor('red')
    
;     xyouts, (wavemax-500)/1e4, min(yr)+(max(yr)-min(yr))*0.8,'z = '+string(zmod[ind_plot[i]],form='(f0.1)'),align=1

  endfor
  xyouts, 0.05,0.5,textoidl('Relative error'),/norm,align=0.5,orien=90

  ps2


;;-- example PCA
  ps1,dir_figs+'egpca.ps'
  i = 0
  !x.margin=[6,3]

  plot, wave_rest_super/1e4,flux_super[*,i]/norm_super[ind_plot[i]], psym=8,xtitle=textoidl('Rest-frame wavelength, /\mum'),ytitle=textoidl('Normalised F_\lambda (with offset)'),$
        xr=[wavemin-500,wavemax+100]/1e4,/xs,yrange=[0,8],/ys,ytickinterval=1
  

  recon = (vwpca_reconstruct(pcs_mock[*,i],evecs[*,0:namplitudes-1])+meanarr)
  oplot, wave_rest_super/1e4,recon,psym=8,color=cgcolor('blue')
  ind = where(flux_mock[*,i] ne 0)
  oplot, ll_eff[ind]/(1+zmod[i])/1e4,flux_mock[ind,i]/norm_super[ind_plot[i]],psym=4
  
  oplot, wave_rest_super/1e4, meanarr+2,psym=8
  oplot, wave_rest_super/1e4, pcs_mock[0,i]*evecs[*,0]+4,psym=8
  oplot, wave_rest_super/1e4, pcs_mock[1,i]*evecs[*,1]+6,psym=8
  oplot, wave_rest_super/1e4, pcs_mock[2,i]*evecs[*,2]+7,psym=8
  xyouts, 1.4, 4.1, 'SC1 = '+string(pcs_mock[0,i],form='(F0.2)'),align=1
  xyouts, 1.4, 6.1, 'SC2 = '+string(pcs_mock[1,i],form='(F0.2)'),align=1
  xyouts, 1.4, 7.1, 'SC3 = '+string(pcs_mock[2,i],form='(F0.2)'),align=1
 

  ps2
 
  splog, 'EGPCA plot model galaxy redshift: ', zmod_phot[i]
;;------------------------------------------------------------------
;;-- investigate model physical parameters
;;------------------------------------------------------------------

  if n_elements(params) ne n_elements(pcs[0,*]) then begin
     print, 'VWSC_FIGS_MODELS: PROBLEM!! pcs do not match params, please fix and rerun'
     stop
 endif

  c_in_aa = 2.99792d18                  ;in [Ang/s]
  AB0 = 2.5*alog10(4*!DPI*3.0857^2*1d5/3.839)

  f_nu = norm_super*(double(wave_norm)^2/c_in_aa) ;in Lsol/Hz
  ABmag = AB0 - 2.5* alog10(f_nu) -48.6

if n_elements(plotparams) eq 0 then range = [-60,240,-30,25] else begin
   pc1range = plotparams.pc1range
   pc2range = plotparams.pc2range
   range = [pc1range[0], pc1range[1], pc2range[0], pc2range[1]]
endelse
if n_elements(plotparams) eq 0 then deltasc = [2.5,1] else deltasc = [plotparams.deltapc1,plotparams.deltapc2]


plot_params,dir_figs+'params_pc12.ps', pcs[0,*],pcs[1,*],params,ABmag,'SC1','SC2',range,deltasc,loadct=13, /ctreverse

if n_elements(plotparams) eq 0 then range=[-60,240,-8,8] else begin
   pc1range = plotparams.pc1range
   pc3range = plotparams.pc3range
   range = [pc1range[0], pc1range[1], pc3range[0],pc3range[1]]
endelse
if n_elements(plotparams) eq 0 then deltasc = [2.5,0.25] else deltasc = [plotparams.deltapc1,plotparams.deltapc3]

plot_params,dir_figs+'params_pc13.ps', pcs[0,*],pcs[2,*],params,ABmag,'SC1','SC3',range,deltasc,loadct=13, /ctreverse

;;-- SED dependence of information
ps1,dir_figs+'deltapcs_agewr.ps'
multiplot,[1,3]
dpc = pcs_mock_zphot[0,*]-pcs[0,*]
plot, alog10(params.age_wr),dpc,psym=3,ytitle=textoidl('\Delta(SC1)'),/xs,yr=[-1.54,1.54],/ys,xr=[7.5,9.7]
hist = histogram(alog10(params.age_wr),min=7.5,max=9.7,binsize=0.1,reverse_ind=R,loc=x) & nbin = n_elements(hist)
per = fltarr(nbin,3)
for i=0,nbin-1 do IF R[i] NE R[i+1] THEN per[i,*] = vw_percentile(dpc[R[R[I] : R[i+1]-1]])
oplot, x+(x[1]-x[0])/2.,per[*,0],color=cgcolor('cyan'),linestyle=2,thick=6
oplot, x+(x[1]-x[0])/2.,per[*,1],color=cgcolor('cyan'),thick=6
oplot, x+(x[1]-x[0])/2.,per[*,2],color=cgcolor('cyan'),linestyle=2,thick=6
multiplot

dpc = pcs_mock_zphot[1,*]-pcs[1,*]
plot, alog10(params.age_wr),dpc,psym=3,ytitle=textoidl('\Delta(SC2)'),/xs,yr=[-1.54,1.54],/ys,xr=[7.5,9.7]
per = fltarr(nbin,3)
for i=0,nbin-1 do IF R[i] NE R[i+1] THEN per[i,*] = vw_percentile(dpc[R[R[I] : R[i+1]-1]])
oplot, x+(x[1]-x[0])/2.,per[*,0],color=cgcolor('cyan'),linestyle=2,thick=6
oplot, x+(x[1]-x[0])/2.,per[*,1],color=cgcolor('cyan'),thick=6
oplot, x+(x[1]-x[0])/2.,per[*,2],color=cgcolor('cyan'),linestyle=2,thick=6
multiplot

dpc = pcs_mock_zphot[2,*]-pcs[2,*] 
plot, alog10(params.age_wr),dpc,psym=3,xtitle=textoidl('log_{10}(age_{wr}/yr)'),ytitle=textoidl('\Delta(SC3)'),/xs,yr=[-1.54,1.54],/ys,xr=[7.5,9.7]
per = fltarr(nbin,3)
for i=0,nbin-1 do IF R[i] NE R[i+1] THEN per[i,*] = vw_percentile(dpc[R[R[I] : R[i+1]-1]])
oplot, x+(x[1]-x[0])/2.,per[*,0],color=cgcolor('cyan'),linestyle=2,thick=6
oplot, x+(x[1]-x[0])/2.,per[*,1],color=cgcolor('cyan'),thick=6
oplot, x+(x[1]-x[0])/2.,per[*,2],color=cgcolor('cyan'),linestyle=2,thick=6
multiplot,/reset
ps2



splog, /close




END
