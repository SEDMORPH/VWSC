FUNCTION ztrend, X,P
RETURN, P[0]*X^P[1]
END



;;****************************************************************************************
;;****************************************************************************************

PRO VWSC_FIGS_DATA, dir_figs, dir_filters, pcs_data, pcerr_data,evecfile,vmaxfile, lines, zphot, minz,maxz, nzbin, lgm, modeltag, nfile, filterset,ind_select, maxmag_survey,mag_obs,area_survey,minmag_survey = minmag_survey

 
  restore, evecfile             ;pcs for models
  nbootstrap = 100
  seed = 234076

;;------------------------------------------------------------------
;;-- Calculate MF and fractions 
;;------------------------------------------------------------------
  if n_elements(minmag_survey) eq 0 then minmag_survey = 14
     splog, 'survey area / deg^2: ',area_survey
     fraction_of_4pi = area_survey/((180/!pi)^2)
     
     LF_nbins=10 & LF_min = -24.5 & LF_max = -20
     MF_nbins=15 & MF_min = 9     & MF_max = 12
     
     complim=0.8

     n_psb_9p75 = (n_psb_10=( n_red_9p75 = (n_red_10=(vsurvey=fltarr(nzbin)))))
     n_psb_9p75_err = (n_psb_10_err= ( n_red_9p75_err = (n_red_10_err=fltarr(nzbin))))
     n_all_9p75 = (n_all_10=( n_sf_9p75 = (n_sf_10=fltarr(nzbin))))
     m_psb_9p75 = (m_psb_10=fltarr(nzbin))
     m_all_9p75 = (m_all_10=fltarr(nzbin))
     m_red_9p75 = (m_red_10=fltarr(nzbin))
     m_sf_9p75 = (m_sf_10=fltarr(nzbin))
     LF_red = (LF_psb = (LF_sf = fltarr(LF_nbins,nzbin)))
     MF_red = (MF_psb = (MF_sf = fltarr(MF_nbins,nzbin)))
     MF_red_err = (MF_psb_err = (MF_sf_err = fltarr(MF_nbins,nzbin)))
  
     ngal = n_elements(zphot)

     VWSC_CLASSIFY, pcs_data, lines,ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty, ind_lomet,ind_junk, class, stack_name, stack_color

     ind_red_all = [ind_red,ind_lomet]
     ind_sf_all = [ind_sf1,ind_sf2,ind_sf3,ind_dusty]


     vmaxarr = fltarr(ngal)
     zbin = fltarr(nzbin+1)
     for i=0,nzbin do zbin[i] = i*(maxz-minz)/float(nzbin)+minz ;end of last bin


     for i=0,nzbin-1 do begin
        splog, '------------------------------------------------------------------'
        splog, 'zrange: ', zbin[i], zbin[i+1]

        ind_z = where(zphot gt zbin[i] and zphot le zbin[i+1])

        ;; vmax for everything in this z bin
        readcol, filterset, filterfiles,form='(A)',/silent ;list of filter files

        readcol, dir_filters+filterfiles[ind_select],lam_fil,fil,/silent
        
        vmaxarr[ind_z] = vwsc_vmax_bc03(zphot[ind_z], mag_obs[ind_z], 0.1, minmag_survey, maxmag_survey, zbin[i], zbin[i+1],$
                                   modeltag,lgm[ind_z].bestmodelid,lam_fil,fil,spec=spec,wave=wave,nfile=nfile)

    endfor

        save, vmaxarr, file=vmaxfile

     for i=0,nzbin-1 do begin

        ind_z = where(zphot gt zbin[i] and zphot le zbin[i+1])
        ind_psb_z = setintersection(ind_z,ind_psb)
        ind_red_all_z = setintersection(ind_z,ind_red_all)
        ind_sf_all_z = setintersection(ind_z,ind_sf_all)

        ; remove objects with no vmax measurement
        ind = where(vmaxarr[ind_z] le 0.0,nn)
        splog, 'nb objects with no vmax measurement',nn
        if ind[0] ne -1 then begin
           ind_z = setdifference(ind_z,ind)
           ind_psb_z = setdifference(ind_psb_z,ind)
           ind_red_all_z =  setdifference(ind_red_all_z,ind)
           ind_sf_all_z =  setdifference(ind_sf_all_z,ind)
        endif

        ;; survey volume
        vsurvey[i] = volume_comoving(zbin[i],zbin[i+1], area=fraction_of_4pi)
        splog, 'survey volume / Mpc^3 : ',vsurvey[i]

        ;; mass functions
        MF_red[*,i] = udspsb_buildLF(lgm[ind_red_all_z].per50, 1/vmaxarr[ind_red_all_z], complim=complim,nbins=MF_nbins,min=MF_min,max=MF_max,loc=MF_x,error=error)-alog10(vsurvey[i])
        MF_red_err[*,i] = error 
        splog, 'Red completeness: ',min(MF_x[where(finite(MF_red[*,i]) eq 1)])
        MF_psb[*,i] = udspsb_buildLF(lgm[ind_psb_z].per50, 1/vmaxarr[ind_psb_z], complim=complim,nbins=MF_nbins,min=MF_min,max=MF_max,error=error)-alog10(vsurvey[i])
        MF_psb_err[*,i] = error 
        splog, 'PSB completeness: ',min(MF_x[where(finite(MF_psb[*,i]) eq 1)])
        MF_sf[*,i] = udspsb_buildLF(lgm[ind_sf_all_z].per50, 1/vmaxarr[ind_sf_all_z], complim=complim,nbins=MF_nbins,min=MF_min,max=MF_max,loc=MF_x,error=error)-alog10(vsurvey[i])
        MF_sf_err[*,i] = error 
        splog, 'SF completeness: ',min(MF_x[where(finite(MF_sf[*,i]) eq 1)])

        ;; number densities, Vmax weighted
        splog, 'nb PSBs  ', n_elements(ind_psb_z)
        splog, 'nb PSBs M* > 10^9.75 ', n_elements(setintersection(ind_psb_z,where(lgm.per50 gt 9.75)))
        splog, 'nb PSBs / Mpc^3 ', n_elements(ind_psb_z)/vsurvey[i]
        splog, 'Volume corrected nb PSBs M*curr > 10^9.75 / Mpc^3 ', total(1/vmaxarr[setintersection(ind_psb_z,where(lgm.per50 gt 9.75))])/vsurvey[i]
        splog, 'Volume corrected nb PSBs M*curr > 10^10 / Mpc^3 ', total(1/vmaxarr[setintersection(ind_psb_z,where(lgm.per50 gt 10))])/vsurvey[i]

        n_all_9p75[i] = total(1/vmaxarr[setintersection(ind_z,where(lgm.per50 gt 9.75))])/vsurvey[i]
        n_all_10[i] = total(1/vmaxarr[setintersection(ind_z,where(lgm.per50 gt 10))])/vsurvey[i]

        n_sf_9p75[i] = total(1/vmaxarr[setintersection(ind_sf_all_z,where(lgm.per50 gt 9.75))])/vsurvey[i]
        n_sf_10[i] = total(1/vmaxarr[setintersection(ind_sf_all_z,where(lgm.per50 gt 10))])/vsurvey[i]
        
        n_psb_9p75[i] = total(1/vmaxarr[setintersection(ind_psb_z,where(lgm.per50 gt 9.75))])/vsurvey[i]
        n_psb_10[i] = total(1/vmaxarr[setintersection(ind_psb_z,where(lgm.per50 gt 10))])/vsurvey[i]
 
        n_red_9p75[i] = total(1/vmaxarr[setintersection(ind_red_all_z,where(lgm.per50 gt 9.75))])/vsurvey[i]
        n_red_10[i] = total(1/vmaxarr[setintersection(ind_red_all_z,where(lgm.per50 gt 10))])/vsurvey[i]
 
        ;; Bootstrap errors on number densities
        ;; this could be done by integrating under the 2D gaussians
        ;; for each point. 
        n_psb_10_bs = (n_red_10_bs = fltarr(nbootstrap))

        for j=0,nbootstrap-1 do begin
           pcs_data_resamp = pcs_data[*,ind_z]
           pcs_data_resamp[0,*] = pcs_data[0,ind_z]+randomn(seed,ngal)*pcerr_data[0,ind_z] ;problem if data=-999 and error=-999 could end up with sensible resamp!
           pcs_data_resamp[1,*] = pcs_data[1,ind_z]+randomn(seed,ngal)*pcerr_data[1,ind_z]
           pcs_data_resamp[2,*] = pcs_data[2,ind_z]+randomn(seed,ngal)*pcerr_data[2,ind_z]
           
           ind = where(finite(pcs_data_resamp) eq 0)
           if ind[0] ne -1 then pcs_data_resamp[ind]=-999
           
           VWSC_CLASSIFY, pcs_data_resamp, lines,ind_red_bs,ind_sf2_bs, ind_sf1_bs, ind_psb_bs, ind_sf3_bs, ind_dusty_bs, ind_lomet_bs,/silent
           ind_red_all_bs = [ind_red_bs,ind_lomet_bs]
           ind_sf_all_bs = [ind_sf1_bs,ind_sf2_bs,ind_sf3_bs,ind_dusty_bs]

           ind = setintersection(ind_psb_bs,where(lgm[ind_z].per50 gt 10))
           n_psb_10_bs[j] = total(1/vmaxarr[ind_z[ind]])/vsurvey[i]

           ind = setintersection(ind_red_all_bs,where(lgm[ind_z].per50 gt 10))
           n_red_10_bs[j] = total(1/vmaxarr[ind_z[ind]])/vsurvey[i]
        endfor
        n_psb_10_err[i] = sqrt((moment(n_psb_10_bs))[1])
        n_red_10_err[i] = sqrt((moment(n_red_10_bs))[1])

        ;; total masses in each class
        ind = setintersection(ind_psb_z,where(lgm.per50 gt 9.75))
        m_psb_9p75[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
        ind = setintersection(ind_psb_z,where(lgm.per50 gt 10))
        m_psb_10[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
        
        ind = setintersection(ind_z,where(lgm.per50 gt 9.75))
        m_all_9p75[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
        ind = setintersection(ind_z,where(lgm.per50 gt 10))
        m_all_10[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
        
        ind = setintersection(ind_red_all_z,where(lgm.per50 gt 9.75))
        m_red_9p75[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
        ind = setintersection(ind_red_all_z,where(lgm.per50 gt 10))
        m_red_10[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
   
        ind = setintersection(ind_sf_all_z,where(lgm.per50 gt 9.75))
        m_sf_9p75[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
        ind = setintersection(ind_sf_all_z,where(lgm.per50 gt 10))
        m_sf_10[i] = total(10^lgm[ind].per50/vmaxarr[ind])/vsurvey[i]
   
         
        ;; luminosity and mass function
;     LF_red[*,i] = udspsb_buildLF(mag_abs_1m[ind_red_all_z], 1/vmaxarr[ind_red_all_z], complim=complim,nbins=LF_nbins,min=LF_min,max=LF_max,loc=LF_x)-alog10(vsurvey[i])
;     LF_psb[*,i] = udspsb_buildLF(mag_abs_1m[ind_psb_z], 1/vmaxarr[ind_psb_z], complim=complim,nbins=LF_nbins,min=LF_min,max=LF_max)-alog10(vsurvey[i])



     endfor

;;-- Calculate fractions
     fraction_psb_9p75 = n_psb_9p75/n_all_9p75
     fraction_red_9p75 = n_red_9p75/n_all_9p75
     fraction_sf_9p75 = n_sf_9p75/n_all_9p75
     fraction_psb_10 = n_psb_10/n_all_10
     fraction_red_10 = n_red_10/n_all_10
     fraction_sf_10 = n_sf_10/n_all_10

     fraction_psb_10_Err = n_psb_10_err/n_all_10 ;assume error on n_all is negligable
     fraction_red_10_Err = n_red_10_err/n_all_10 ;assume error on n_all is negligable


;;-- plot fraction as a function of z
     zz = (zbin[1]-zbin[0])/2.+zbin[0:nzbin-1] ;-1 is correct
     ps1, dir_figs+'fraction_z.ps'
     !x.charsize=(!y.charsize= (!p.charsize=1.2))
     plot,zz,fraction_red_10,/nodata,ytitle='Fraction of quiescent and PSB',xr=minmax(zbin),xtitle='z',/xs
     oploterror,zz,fraction_red_10,fraction_red_10_err,psym=-2,color=cgcolor('red')
     oplot,zz,fraction_red_9p75,psym=-1,color=cgcolor('red')
     
     oploterror,zz,fraction_psb_10,fraction_psb_10_err,psym=-2,color=cgcolor('orange')
     oplot,zz,fraction_psb_9p75,psym=-1,color=cgcolor('orange')
     
     zfit = mpfitfun('ztrend',zz,fraction_psb_10,fraction_psb_10_err,[0.1,3],yfit=yfit,perror=perror)
     oplot,zz,yfit,linestyle=1
     splog, 'Ztrend results: ', zfit
     splog, 'Errors on fit: ', perror
     
     xyouts, 0.7, 0.78, '+ log(M*)>9.75',/norm,align=0
     xyouts, 0.7, 0.85, '* log(M*)>10',/norm,align=0
     ps2

;;-- plot number growth
     ndot_psb_10 = fraction_psb_10/0.5 ;t_PSB = 0.5Gyr ? about right according to topcat models
     ndot_psb_10_err =  fraction_psb_10_err/0.5
     ndot_psb_9p75 = fraction_psb_9p75/0.5
     t_uni = galage(zz,1000)/1d9

     dt_reddot = (t_uni[indgen(nzbin-1)]-t_uni[indgen(nzbin-1)+1])
     ndot_red_10 = ( fraction_red_10[indgen(nzbin-1)]-fraction_red_10[indgen(nzbin-1)+1])/dt_reddot ;ok
     ndot_red_10_err = sqrt(fraction_red_10_err[indgen(nzbin-1)]^2 + (fraction_red_10_err[indgen(nzbin-1)+1])^2)/dt_reddot
     ndot_red_9p75 = ( fraction_red_9p75[indgen(nzbin-1)]-fraction_red_9p75[indgen(nzbin-1)+1])/dt_reddot
;     linterp, zz,n_all_10,zbin[indgen(nzbin-1)+1],n_all_interp_10
;     linterp, zz,n_all_9p75,zbin[indgen(nzbin-1)+1],n_all_interp_9p75

     t_reddot = t_uni[indgen(nzbin-1)+1]+dt_reddot/2.

;; need to mulitply by tot number dens, average over full survey to
;; avoid cosmic variance

     ps1, dir_figs+'number_growth.ps'
     !x.charsize=(!y.charsize= (!p.charsize=1.2))
     plot,zbin[indgen(nzbin-1)+1],ndot_red_10,/nodata,xtitle='z',ytitle=textoidl('dn/dt [Mpc^{-3} Gyr^{-1}]'),title=textoidl('t_{PSB} = 0.5 Gyr'),xr=[minz,maxz],yr=[0,0.2]
     oploterror,zbin[indgen(nzbin-1)+1],ndot_red_10,ndot_red_10_err, psym=-2,color=cgcolor('red')
     oplot,zbin[indgen(nzbin-1)+1],ndot_red_9p75, psym=-1,color=cgcolor('red')

     oploterror,zz,ndot_psb_10,ndot_psb_10_err,psym=-2,color=cgcolor('orange')
     oplot,zz,ndot_psb_9p75,psym=-1,color=cgcolor('orange')

     xyouts, 0.7, 0.78, '+ log(M*)>9.75',/norm,align=0
     xyouts, 0.7, 0.85, '* log(M*)>10',/norm,align=0

     ps2

;;-- plot mass functions
     ps1, dir_figs+'MF.ps'
     device,ysize=9,xsize=25
     multiplot, [nzbin,1],mytitle=textoidl('log \phi / Mpc^{-3} dex^{-1}'),mxtitle=textoidl('Log(M*/M_\odot)')
     nbin = n_elements(mf_x)
     for i=0,nzbin-1 do begin
        plot, MF_x,MF_red[*,i],xr=[9,12],yr=[-5.5,-2.5],/nodata,/ys,/xs
        oploterror, MF_x,MF_red[*,i],MF_red_err[*,i],psym=-1,color=cgcolor('red'),errcol=cgcolor('red')
        oploterror, MF_x,MF_sf[*,i],MF_sf_err[*,i],psym=-1,color=cgcolor('blue'),errcol=cgcolor('blue')
        oploterror, MF_x,MF_psb[*,i],MF_psb_err[*,i],psym=-1,color=cgcolor('orange'),errcol=cgcolor('orange')
        xyouts, 9.7,-2.7,string(zbin[i],form='(F0.1)')+' < z < ' + string(zbin[i+1],form='(F0.1)')
        multiplot
     endfor
     multiplot,/reset
     ps2

;;-- plot normalised mass function
     ps1, dir_figs+'MF_norm.ps'
     device,ysize=9,xsize=25
     multiplot, [nzbin,1],mytitle=textoidl('normalised log \phi / Mpc^{-3} dex^{-1}'),mxtitle=textoidl('Log(M*/M_\odot)')
     nbin = n_elements(mf_x)
     for i=0,nzbin-1 do begin
    ;    ind = where(mf_x gt 9.75)
        plot, MF_x[*],MF_red[*,i]-max(MF_red[*,i]),xr=[9,12],yr=[-2.,0.5],/nodata,/ys,/xs
        oplot, MF_x[*],MF_red[*,i]-max(MF_red[*,i]),psym=-1,color=cgcolor('red')
;        oplot, MF_x[*],MF_sf[*,i]-max(MF_sf[*,i]),psym=-1,color=cgcolor('blue')
        oplot, MF_x[*],MF_psb[*,i]-max(MF_psb[*,i]),psym=-1,color=cgcolor('orange')
        xyouts, 9.7,-1.85,string(zbin[i],form='(F0.1)')+' < z < ' + string(zbin[i+1],form='(F0.1)')
        multiplot
     endfor
     multiplot,/reset
     ps2


END
