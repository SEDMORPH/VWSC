PRO VWSC_PLOT_PARAMS, psfile, modelfile,xtitle, ytitle, range,bins,_extra=extra

  restore, modelfile

;;-- SC1 vs. SC2 or SC1 vs. SC3
  CASE xtitle of
     'SC1': x = pcs[0,*]
     'SC2': x = pcs[1,*]
     'SC3': x = pcs[2,*]
     ELSE : BEGIN 
        print, 'xtitle must be string SC1, SC2 or SC3'
        return
        END
  ENDCASE
  CASE ytitle of
     'SC1': y = pcs[0,*]
     'SC2': y = pcs[1,*]
     'SC3': y = pcs[2,*]
     ELSE : BEGIN 
        print, 'ytitle must be string SC1, SC2 or SC3'
        return
     END
  ENDCASE

  loadct = 13                   ;rainbow

  ps1,psfile
  ncolors=100
  loadct,10
  plotsym,0,0.1,/fill
  !x.margin=[8,3]
  multiplot,[2,2]
  
  nmod = n_elements(params)
  
  vwplot_2dhist,x,y,range,bins, weight1=alog10(params.age_wr),vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out, ytitle=ytitle,/ctreverse
  cgcolorbar,position= [0.46,0.55,0.48,0.90],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
  xyouts, 0.47,0.92,textoidl('log(Age_{wr})'),/norm,align=0.5,charsize=0.8
  
  multiplot
  vwplot_2dhist,x,y,range,bins, weight1=params.zmet,vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out,/ctreverse
  cgcolorbar,position= [0.89,0.55,0.91,0.90],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
  xyouts, 0.93,0.92,textoidl('Z'),/norm,align=0.5,charsize=0.8
  
  multiplot
  vwplot_2dhist,x,y,range,bins, weight1=params.tauv0,vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out, ytitle=ytitle,xtitle=xtitle,/ctreverse

  cgcolorbar,position= [0.46,0.13,0.48,0.48],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.8,divisions=5,/right
  xyouts, 0.49,0.5,textoidl('\tau_v'),/norm,align=0.5,charsize=0.8
  
  multiplot
;; time: 1.e6, 1.e7, 1.e8, 1.e9 and 2.e9
  vwplot_2dhist,x,y,range,bins, weight1=params.fburst[3],vmax = fltarr(nmod)+1,/nobar,loadct=loadct,range_out=range_out,xtitle=xtitle
  cgcolorbar,position= [0.89,0.13,0.91,0.48],range=range_out,/noerase,format='(F0.1)',/vertical,charsize=0.82,divisions=5,/reverse,/right
  xyouts, 0.93,0.5,textoidl('F_{burst}'),/norm,align=0.5,charsize=0.8
  multiplot,/reset
  ps2

END
