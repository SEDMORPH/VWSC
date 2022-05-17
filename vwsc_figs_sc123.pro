PRO fig_pc123, pcs,stack_color, pcs_data, pcerr_data, ind_red, ind_sf1, ind_sf2, ind_sf3,ind_psb,ind_lomet, ind_dusty, $                  ;REQUIRED
               pcs_spec,ind_red_spec, ind_sf1_spec, ind_sf2_spec, ind_sf3_spec,ind_psb_spec,ind_lomet_spec, ind_dusty_spec, $ ;OPTIONAL
               x_psb1=x_psb,y_psb1=y_psb, x_psb2=x_psb2,y_psb2=y_psb2,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty, plotparams=plotparams


  if n_elements(plotparams) eq 0 then begin
     pc1range = [-50,150]
     pc2range = [-50,20]
     pc3range = [-10,15]
     deltapc1 = 2.5
     deltapc2 = 0.5
     deltapc3 = 0.25
  endif else begin
     pc1range = plotparams.pc1range
     pc2range = plotparams.pc2range
     pc3range = plotparams.pc3range
     deltapc1 = plotparams.deltapc1
     deltapc2 = plotparams.deltapc2
     deltapc3 = plotparams.deltapc3
  endelse

;;------------------------------------------------------------------
;;-- Supercolour histogram
;;------------------------------------------------------------------
;  range_in = [5,100]
  ;levels=[0.75,1.2,1.4,1.6,1.75]
  vwplot_2dhist, pcs_data[0,*], pcs_data[1,*],[pc1range[0],pc1range[1],pc2range[0],pc2range[1]],[deltapc1,deltapc2],psym=8,xr=[pc1range[0],pc1range[1]],yr=[pc2range[0],pc2range[1]],/xs,/ys,ytitle='SC2',/log,/nobar,range_in=range_in,loadct=10,range_out=range_out,range_color=range_color
  position = [0.32,0.55,0.33,0.9]
  cgcolorbar, range =alog10([range_out[0],range_out[1]]), /vertical, /right,/reverse,format = '(F0.1)', position=position, charsize = 0.5*!p.charsize,tickinterval=0.2,bottom=min(range_color)
  xyouts, position[0]+0.02, position[3]+0.02, 'log(N)',/normal,align=0.5

  ;; oplot, x_psb2, y_psb2,color=cgcolor('blue')
  ;; oplot, x_psb, y_psb,color=cgcolor('blue')
  ;; oplot, x_red, y_red,color=cgcolor('blue')
  ;; oplot, x_dusty, y_dusty,color=cgcolor('blue')
  ;; oplot, [-10,-10],[-15,3],color=cgcolor('blue')

 

  ;;-- PC12 points + spectra points
  multiplot
  plotsym,0,0.1,/fill
  plot,  pcs_data[0,*],pcs_data[1,*],psym=8,xr=[pc1range[0],pc1range[1]],yr=[pc2range[0],pc2range[1]],/xs,/ys,/nodata
  if n_elements(ind_red) ge 3 then  oplot,  pcs_data[0,ind_red],pcs_data[1,ind_red],color=cgcolor(stack_color[0]),psym=8
  if n_elements(ind_sf2) ge 3 then  oplot,  pcs_data[0,ind_sf2],pcs_data[1,ind_sf2],color=cgcolor(stack_color[2]),psym=8
  if n_elements(ind_sf3) ge 3 then   oplot,  pcs_data[0,ind_sf3],pcs_data[1,ind_sf3],color=cgcolor(stack_color[1]),psym=8
  if n_elements(ind_sf1) ge 3 then   oplot,  pcs_data[0,ind_sf1],pcs_data[1,ind_sf1],color=cgcolor(stack_color[3]),psym=8
  if n_elements(ind_dusty) ge 3 then   oplot,  pcs_data[0,ind_dusty],pcs_data[1,ind_dusty],color=cgcolor(stack_color[5]),psym=8

;  plotsym,0,0.3 & oplot,  pcs_data[0,ind_psb],pcs_data[1,ind_psb],color=cgcolor('black'),psym=8
  plotsym,0,0.1,/fill & if ind_psb[0] ne -1 then oplot,  pcs_data[0,ind_psb],pcs_data[1,ind_psb],color=cgcolor(stack_color[4]),psym=8
;  plotsym,0,0.3 & if n_elements(ind_lomet) ne 0 then oplot,  pcs_data[0,ind_lomet],pcs_data[1,ind_lomet],color=cgcolor('black'),psym=8
  plotsym,0,0.1,/fill & if ind_lomet[0] ne -1 then oplot,  pcs_data[0,ind_lomet],pcs_data[1,ind_lomet],color=cgcolor(stack_color[6]),psym=8

  ;; labels
  plotsym,0,0.1,/fill & plots, 0.58,0.94,psym=8,/norm,color=cgcolor(stack_color[3]) & xyouts, 0.59, 0.93,'SF1',/norm,color=cgcolor(stack_color[3])
  plotsym,0,0.1,/fill & plots, 0.58,0.94-0.02,psym=8,/norm,color=cgcolor(stack_color[2]) & xyouts, 0.59, 0.93-0.02,'SF2',/norm,color=cgcolor(stack_color[2])
  plotsym,0,0.1,/fill & plots, 0.58,0.94-0.04,psym=8,/norm,color=cgcolor(stack_color[1]) & xyouts, 0.59, 0.93-0.04,'SF3',/norm,color=cgcolor(stack_color[1])
  plotsym,0,0.1,/fill & plots, 0.58,0.94-0.06,psym=8,/norm,color=cgcolor(stack_color[0]) & xyouts, 0.59, 0.93-0.06,'Red',/norm,color=cgcolor(stack_color[0])
  plotsym,0,0.3 & plots, 0.58,0.94-0.08,psym=8,/norm,color=cgcolor('black') & plotsym,0,0.2,/fill & plots, 0.58,0.94-0.08,psym=8,/norm,color=cgcolor(stack_color[4]) & xyouts, 0.59, 0.93-0.08,'PSB',/norm,color=cgcolor(stack_color[4])
  plotsym,0,0.1,/fill & plots, 0.58,0.94-0.1,psym=8,/norm,color=cgcolor(stack_color[5]) & xyouts, 0.59, 0.93-0.1,'Dusty',/norm,color=cgcolor(stack_color[5])
  plotsym,0,0.3 & plots, 0.58,0.94-0.12,psym=8,/norm,color=cgcolor('black') & plotsym,0,0.2,/fill & plots, 0.58,0.94-0.12,psym=8,/norm,color=cgcolor(stack_color[6]) & xyouts, 0.59, 0.93-0.12,'Low-met',/norm,color=cgcolor(stack_color[6])

  ;; typical errors
  if n_elements(ind_red) ge 3 then  oploterror, -40, -40,median(pcerr_data[0,ind_red]),median(pcerr_data[1,ind_red]),errcolor=cgcolor(stack_color[0])
  if n_elements(ind_lomet) ge 3 then oploterror, -20, -40,median(pcerr_data[0,ind_lomet]),median(pcerr_data[1,ind_lomet]),errcolor=cgcolor(stack_color[6])
  if n_elements(ind_dusty) ge 3 then  oploterror, 0, -40,median(pcerr_data[0,ind_dusty]),median(pcerr_data[1,ind_dusty]),errcolor=cgcolor(stack_color[5])
  if n_elements(ind_psb) ge 3 then  oploterror, 20, -40,median(pcerr_data[0,ind_psb]),median(pcerr_data[1,ind_psb]),errcolor=cgcolor(stack_color[4])
  if n_elements(ind_sf3) ge 3 then  oploterror, 40, -40,median(pcerr_data[0,ind_sf3]),median(pcerr_data[1,ind_sf3]),errcolor=cgcolor(stack_color[1])
  if n_elements(ind_sf2) ge 3 then  oploterror, 60, -40,median(pcerr_data[0,ind_sf2]),median(pcerr_data[1,ind_sf2]),errcolor=cgcolor(stack_color[2])
  if n_elements(ind_sf1) ge 3 then  oploterror, 80, -40,median(pcerr_data[0,ind_sf1]),median(pcerr_data[1,ind_sf1]),errcolor=cgcolor(stack_color[3])

;  plotsym,0,0.6,/fill
;  oplot,  pcs_spec[0,ind_red_spec],pcs_spec[1,ind_red_spec],color=cgcolor(stack_color[0]),psym=8
;  oplot,  pcs_spec[0,ind_sf2_spec],pcs_spec[1,ind_sf2_spec],color=cgcolor(stack_color[2]),psym=8
;  oplot,  pcs_spec[0,ind_sf3_spec],pcs_spec[1,ind_sf3_spec],color=cgcolor(stack_color[1]),psym=8
;  oplot,  pcs_spec[0,ind_sf1_spec],pcs_spec[1,ind_sf1_spec],color=cgcolor(stack_color[3]),psym=8
;  oplot,  pcs_spec[0,ind_dusty_spec],pcs_spec[1,ind_dusty_spec],color=cgcolor(stack_color[5]),psym=8
;  if ind_psb[0] ne -1 then oplot,  pcs_spec[0,ind_psb_spec],pcs_spec[1,ind_psb_spec],color=cgcolor(stack_color[4]),psym=8
 

  ;;-- PC12 models + spectra points 
  multiplot
  vwplot_2dhist, pcs[0,*],pcs[1,*],[pc1range[0],pc1range[1],pc2range[0],pc2range[1]],[deltapc1,deltapc2],psym=8,xr=[pc1range[0],pc1range[1]],yr=[pc2range[0],pc2range[1]],/xs,/ys,/log,/nobar,loadct=0,range_out=range_out,range_color=range_color,/nodata
  position = [0.92,0.55,0.93,0.9]
  cgcolorbar, range =alog10([range_out[0],range_out[1]]), /vertical, /right,/reverse,format = '(F0.1)', position=position, charsize = 0.5*!p.charsize,tickinterval=0.5
  xyouts, position[0]+0.02, position[3]+0.02, 'log(N)',/normal,align=0.5

  if n_elements(pcs_spec) ne 0 then begin

     plotsym,0,0.1,/fill
     if ind_red_spec[0] ne -1 then oplot,  pcs_spec[0,ind_red_spec],pcs_spec[1,ind_red_spec],color=cgcolor(stack_color[0]),psym=8
     if ind_sf2_spec[0] ne -1 then oplot,  pcs_spec[0,ind_sf2_spec],pcs_spec[1,ind_sf2_spec],color=cgcolor(stack_color[2]),psym=8
     if ind_sf3_spec[0] ne -1 then oplot,  pcs_spec[0,ind_sf3_spec],pcs_spec[1,ind_sf3_spec],color=cgcolor(stack_color[1]),psym=8
     if ind_sf1_spec[0] ne -1 then  oplot,  pcs_spec[0,ind_sf1_spec],pcs_spec[1,ind_sf1_spec],color=cgcolor(stack_color[3]),psym=8
     if ind_dusty_spec[0] ne -1 then oplot,  pcs_spec[0,ind_dusty_spec],pcs_spec[1,ind_dusty_spec],color=cgcolor(stack_color[5]),psym=8
     if ind_psb_spec[0] ne -1 then oplot,  pcs_spec[0,ind_psb_spec],pcs_spec[1,ind_psb_spec],color=cgcolor(stack_color[4]),psym=8
     if ind_lomet_spec[0] ne -1 then oplot,  pcs_spec[0,ind_lomet_spec],pcs_spec[1,ind_lomet_spec],color=cgcolor(stack_color[6]),psym=8

  endif


  ;;-- PC13 histogram
  multiplot
;  range_in = [5,100]
;  levels=[0.75,1.25,1.5,1.75],
  vwplot_2dhist, pcs_data[0,*],pcs_data[2,*],[pc1range[0],pc1range[1],pc3range[0],pc3range[1]],[deltapc1,deltapc3],psym=8,xr=[pc1range[0],pc1range[1]],yr=[pc3range[0],pc3range[1]],/xs,/ys,ytitle='SC3',xtitle='SC1',/log,/nobar,range_in=range_in,loadct=10,range_out=range_out,range_color=range_color
  position = [0.32,0.11,0.33,0.45]
  cgcolorbar, range =alog10([range_out[0],range_out[1]]), /vertical, /right,/reverse,format = '(F0.1)', position=position, charsize = 0.5*!p.charsize,tickinterval=0.2,bottom=min(range_color)
  xyouts, position[0]+0.02, position[3]+0.02, 'log(N)',/normal,align=0.5

;  oplot, [-10,-10],[-15,3],color=cgcolor('blue')


  ;;-- PC13 points + spectra points
  multiplot
  plotsym,0,0.1,/fill
  plot,  pcs_data[0,*],pcs_data[2,*],psym=3,xr=[pc1range[0],pc1range[1]],yr=[pc3range[0],pc3range[1]],/xs,/ys,xtitle='SC1',/nodata
   if n_elements(ind_red) ge 3 then  oplot,  pcs_data[0,ind_red],pcs_data[2,ind_red],color=cgcolor(stack_color[0]),psym=8
   if n_elements(ind_sf2) ge 3 then  oplot,  pcs_data[0,ind_sf2],pcs_data[2,ind_sf2],color=cgcolor(stack_color[2]),psym=8
   if n_elements(ind_sf3) ge 3 then  oplot,  pcs_data[0,ind_sf3],pcs_data[2,ind_sf3],color=cgcolor(stack_color[1]),psym=8
   if n_elements(ind_sf1) ge 3 then  oplot,  pcs_data[0,ind_sf1],pcs_data[2,ind_sf1],color=cgcolor(stack_color[3]),psym=8
  if n_elements(ind_dusty) ge 3 then   oplot,  pcs_data[0,ind_dusty],pcs_data[2,ind_dusty],color=cgcolor(stack_color[5]),psym=8
  
  
;  plotsym,0,0.3 & oplot,  pcs_data[0,ind_psb],pcs_data[2,ind_psb],color=cgcolor('black'),psym=8
  plotsym,0,0.1,/fill & if ind_psb[0] ne -1 then oplot,  pcs_data[0,ind_psb],pcs_data[2,ind_psb],color=cgcolor(stack_color[4]),psym=8
;  plotsym,0,0.3 &  if n_elements(ind_lomet) ne 0 then oplot,  pcs_data[0,ind_lomet],pcs_data[2,ind_lomet],color=cgcolor('black'),psym=8
  plotsym,0,0.1,/fill &  if ind_lomet[0] ne -1 then oplot,  pcs_data[0,ind_lomet],pcs_data[2,ind_lomet],color=cgcolor(stack_color[6]),psym=8

   ;; typical errors
  if n_elements(ind_red) ge 3 then  oploterror, -40, 10,median(pcerr_data[0,ind_red]),median(pcerr_data[2,ind_red]),errcolor=cgcolor(stack_color[0])
  if n_elements(ind_lomet) ge 3 then  oploterror, -20, 10,median(pcerr_data[0,ind_lomet]),median(pcerr_data[2,ind_lomet]),errcolor=cgcolor(stack_color[6])
  if n_elements(ind_dusty) ge 3 then  oploterror, 0, 10,median(pcerr_data[0,ind_dusty]),median(pcerr_data[2,ind_dusty]),errcolor=cgcolor(stack_color[5])
  if n_elements(ind_psb) ge 3 then  oploterror, 20, 10,median(pcerr_data[0,ind_psb]),median(pcerr_data[2,ind_psb]),errcolor=cgcolor(stack_color[4])
  if n_elements(ind_sf3) ge 3 then  oploterror, 40, 10,median(pcerr_data[0,ind_sf3]),median(pcerr_data[2,ind_sf3]),errcolor=cgcolor(stack_color[1])
  if n_elements(ind_sf2) ge 3 then  oploterror, 60, 10,median(pcerr_data[0,ind_sf2]),median(pcerr_data[2,ind_sf2]),errcolor=cgcolor(stack_color[2])
  if n_elements(ind_sf1) ge 3 then  oploterror, 80, 10,median(pcerr_data[0,ind_sf1]),median(pcerr_data[2,ind_sf1]),errcolor=cgcolor(stack_color[3])


;  plotsym,0,0.6,/fill
;  oplot,  pcs_spec[0,ind_red],pcs_spec[2,ind_red],color=cgcolor(stack_color[0]),psym=8
;  oplot,  pcs_spec[0,ind_sf2],pcs_spec[2,ind_sf2],color=cgcolor(stack_color[2]),psym=8
;  oplot,  pcs_spec[0,ind_sf3],pcs_spec[2,ind_sf3],color=cgcolor(stack_color[1]),psym=8
;  oplot,  pcs_spec[0,ind_sf1],pcs_spec[2,ind_sf1],color=cgcolor(stack_color[3]),psym=8
;  oplot,  pcs_spec[0,ind_dusty],pcs_spec[2,ind_dusty],color=cgcolor(stack_color[5]),psym=8
;  if ind_psb[0] ne -1 then oplot,  pcs_spec[0,ind_psb],pcs_spec[2,ind_psb],color=cgcolor(stack_color[4]),psym=8


  ;;-- PC13 models + spectra points
  multiplot
  vwplot_2dhist, pcs[0,*],pcs[2,*],[pc1range[0],pc1range[1],pc3range[0],pc3range[1]],[deltapc1,deltapc3],psym=8,xr=[pc1range[0],pc1range[1]],yr=[pc3range[0],pc3range[1]],/xs,/ys,xtitle='SC1',/log,/nobar,loadct=0,range_out=range_out,range_color=range_color,/nodata
  position = [0.92,0.11,0.93,0.45]
  cgcolorbar, range =alog10([range_out[0],range_out[1]]), /vertical, /right,/reverse,format = '(F0.1)', position=position, charsize = 0.5*!p.charsize,tickinterval=0.5
  xyouts, position[0]+0.02, position[3]+0.02, 'log(N)',/normal,align=0.5

  if n_elements(pcs_spec) ne 0 then begin
     plotsym,0,0.1,/fill
     if ind_red_spec[0] ne -1 then oplot,  pcs_spec[0,ind_red_spec],pcs_spec[2,ind_red_spec],color=cgcolor(stack_color[0]),psym=8
     if ind_sf2_spec[0] ne -1 then oplot,  pcs_spec[0,ind_sf2_spec],pcs_spec[2,ind_sf2_spec],color=cgcolor(stack_color[2]),psym=8
     if ind_sf3_spec[0] ne -1 then oplot,  pcs_spec[0,ind_sf3_spec],pcs_spec[2,ind_sf3_spec],color=cgcolor(stack_color[1]),psym=8
     if ind_sf1_spec[0] ne -1 then oplot,  pcs_spec[0,ind_sf1_spec],pcs_spec[2,ind_sf1_spec],color=cgcolor(stack_color[3]),psym=8
     if ind_dusty_spec[0] ne -1 then oplot,  pcs_spec[0,ind_dusty_spec],pcs_spec[2,ind_dusty_spec],color=cgcolor(stack_color[5]),psym=8
     if ind_psb_spec[0] ne -1 then oplot,  pcs_spec[0,ind_psb_spec],pcs_spec[2,ind_psb_spec],color=cgcolor(stack_color[4]),psym=8
     if ind_lomet_spec[0] ne -1 then oplot,  pcs_spec[0,ind_lomet_spec],pcs_spec[2,ind_lomet_spec],color=cgcolor(stack_color[6]),psym=8
  endif


end

;;****************************************************************************************
;;****************************************************************************************

PRO VWSC_FIGS_SC123, dir_figs, dir_filters, pcs_data,pcerr_data, evecfile, lines, zphot, minz,maxz, nzbin, pcs_spec=pcs_spec, zspec=zspec,plotparams=plotparams


  restore, evecfile             ;pcs for models


  ps1, dir_figs+'PC123_class.ps'
  device,/portrait,xoffset=0,yoffset=0,ysize=20,xsize=30
  multiplot, [3,2],ygap=0.01,xgap=0.01,mytitoffset=-1,/doyaxis

  VWSC_CLASSIFY, pcs_data, lines, ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty, ind_lomet, ind_junk,class, stack_name, stack_color, x_psb1=x_psb,y_psb1=y_psb, x_psb2=x_psb2,y_psb2=y_psb2,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty
  if n_elements(pcs_spec) ne 0 then VWSC_CLASSIFY, pcs_spec, lines, ind_red_spec,ind_sf2_spec, ind_sf1_spec, ind_psb_spec, ind_sf3_spec, ind_dusty_spec, ind_lomet_spec

  fig_pc123, pcs, stack_color, pcs_data, pcerr_data,ind_red, ind_sf1, ind_sf2, ind_sf3, ind_psb, ind_lomet, ind_dusty, pcs_spec, ind_red_spec, ind_sf1_spec, ind_sf2_spec, ind_sf3_spec,ind_psb_spec,ind_lomet_spec, ind_dusty_spec, x_psb1=x_psb,y_psb1=y_psb, x_psb2=x_psb2,y_psb2=y_psb2,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty,plotparams=plotparams
  

  multiplot,/reset
  ps2



;;--- PC123 plots as function z

  if n_elements(nzbin) ne 0 then begin
     ps1, dir_figs+'PC123_zbin.ps'
     device,/portrait,xoffset=0,yoffset=0,ysize=20,xsize=30

     zbin= fltarr(nzbin+1)
     for i=0,nzbin-1 do begin
        multiplot, [3,2],ygap=0.01,xgap=0.01,mytitoffset=-1,/doyaxis,mtitle=string(i*(maxz-minz)/float(nzbin)+minz,form='(F0.2)')+'< z <'+string((i+1)*(maxz-minz)/float(nzbin)+minz,form='(F0.2)')
        
        ind = where(zphot gt i*(maxz-minz)/float(nzbin)+minz and zphot le  (i+1)*(maxz-minz)/float(nzbin)+minz)
        zbin[i] = i*(maxz-minz)/float(nzbin)+minz

        VWSC_CLASSIFY, pcs_data[*,ind], lines,ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty, ind_lomet,ind_junk, class, stack_name, stack_color
          
        if n_elements(pcs_spec) ne 0 then begin
           ind_spec = where(zspec gt i*(maxz-minz)/float(nzbin)+minz and zspec le  (i+1)*(maxz-minz)/float(nzbin)+minz)
           VWSC_CLASSIFY, pcs_spec[*,ind_spec], lines, ind_red_spec,ind_sf2_spec, ind_sf1_spec, ind_psb_spec, ind_sf3_spec, ind_dusty_spec, ind_lomet_spec

           fig_pc123,  pcs, stack_color,pcs_data[*,ind],pcerr_data[*,ind], ind_red, ind_sf1, ind_sf2, ind_sf3,ind_psb,ind_lomet, ind_dusty,pcs_spec[*,ind_spec],ind_red_spec, ind_sf1_spec, ind_sf2_spec, ind_sf3_spec,ind_psb_spec,ind_lomet_spec, ind_dusty_spec, x_psb1=x_psb,y_psb1=y_psb, x_psb2=x_psb2,y_psb2=y_psb2,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty,plotparams=plotparams
        endif else begin
           fig_pc123,  pcs, stack_color,pcs_data[*,ind],pcerr_data[*,ind], ind_red, ind_sf1, ind_sf2, ind_sf3,ind_psb,ind_lomet, ind_dusty, x_psb1=x_psb,y_psb1=y_psb, x_psb2=x_psb2,y_psb2=y_psb2,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty,plotparams=plotparams
        endelse

        
        
        erase
     endfor
     zbin[i] = i*(maxz-minz)/float(nzbin)+minz ;end of last bin

     ps2
 
;;--- PC123 plots as function time

     time = galage(zphot,1000)/1e9
     mint = min(time)
     maxt = max(time)

     ps1, dir_figs+'PC123_tbin.ps'
     device,/portrait,xoffset=0,yoffset=0,ysize=20,xsize=30

     tt = fltarr(nzbin)
     for i=0,nzbin-1 do begin
        multiplot, [3,2],ygap=0.01,xgap=0.01,mytitoffset=-1,/doyaxis,mtitle=string(i*(maxt-mint)/float(nzbin)+mint,form='(F0.2)')+'< t <'+string((i+1)*(maxt-mint)/float(nzbin)+mint,form='(F0.2)')
        
        ind = where(time gt i*(maxt-mint)/float(nzbin)+mint and time le  (i+1)*(maxt-mint)/float(nzbin)+mint)
        tt[i] = ((i+1)*(maxt-mint)/float(nzbin) -i*(maxt-mint)/float(nzbin))/2.
     
        VWSC_CLASSIFY, pcs_data[*,ind], lines,ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty, ind_lomet,ind_junk, class, stack_name, stack_color
        
        if n_elements(pcs_spec) ne 0 then begin
           ind_spec = where(zspec gt i*(maxz-minz)/float(nzbin)+minz and zspec le  (i+1)*(maxz-minz)/float(nzbin)+minz)
           VWSC_CLASSIFY, pcs_spec[*,ind_spec], lines, ind_red_spec,ind_sf2_spec, ind_sf1_spec, ind_psb_spec, ind_sf3_spec, ind_dusty_spec, ind_lomet_spec
           fig_pc123,  pcs, stack_color,pcs_data[*,ind],pcerr_data[*,ind], ind_red, ind_sf1, ind_sf2, ind_sf3,ind_psb,ind_lomet, ind_dusty,pcs_spec[*,ind_spec],ind_red_spec, ind_sf1_spec, ind_sf2_spec, ind_sf3_spec,ind_psb_spec,ind_lomet_spec, ind_dusty_spec, x_psb1=x_psb,y_psb1=y_psb, x_psb2=x_psb2,y_psb2=y_psb2,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty,plotparams=plotparams
        endif else fig_pc123,  pcs, stack_color,pcs_data[*,ind],pcerr_data[*,ind], ind_red, ind_sf1, ind_sf2, ind_sf3,ind_psb,ind_lomet, ind_dusty, x_psb1=x_psb,y_psb1=y_psb, x_psb2=x_psb2,y_psb2=y_psb2,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty,plotparams=plotparams

        
        
        erase
     endfor
     
     ps2
  
  endif


END
