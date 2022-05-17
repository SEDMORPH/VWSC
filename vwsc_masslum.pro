PRO VWSC_MASSLUM




;;------------------------------------------------------------------
;;-- Mass, colour, luminosity, UVJ
;;------------------------------------------------------------------
  if n_elements(massfile) ne 0 then begin
     ps1,dir_figs+'mag_colour_data.ps'
     !p.multi=0
     plotsym,0,0.2,/fill
     plot, mag_abs_1m, uvj[0,*]-uvj[1,*],psym=3,ytitle='U-V',xtitle=textoidl('M_{AB,1\mum}'),xr=[-20,-24],yr=[0.1,2.5],/nodata,/ys
     oplot, mag_abs_1m[ind_red],uvj[0,ind_red]-uvj[1,ind_red],psym=8,color=cgcolor(stack_color[0])
     oplot,mag_abs_1m[ind_sf2], uvj[0,ind_sf2]-uvj[1,ind_sf2],psym=8,color=cgcolor(stack_color[2])
     oplot, mag_abs_1m[ind_sb],uvj[0,ind_sb]-uvj[1,ind_sb],psym=8,color=cgcolor(stack_color[3])
     oplot, mag_abs_1m[ind_sf3],uvj[0,ind_sf3]-uvj[1,ind_sf3],psym=8,color=cgcolor(stack_color[1])
;  plotsym,0,0.4,/fill
     oplot, mag_abs_1m[ind_dusty],uvj[0,ind_dusty]-uvj[1,ind_dusty],psym=8,color=cgcolor(stack_color[5])
     oplot,mag_abs_1m[ind_psb], uvj[0,ind_psb]-uvj[1,ind_psb],psym=8,color=cgcolor(stack_color[4])
     oplot,mag_abs_1m[ind_lomet], uvj[0,ind_lomet]-uvj[1,ind_lomet],psym=8,color=cgcolor(stack_color[6])

     x = findgen(10)-23
     y = -0.1*x-0.7
     oplot, x,y,linestyle=2

     y =-0.13*x-1.7
     oplot, x,y,linestyle=2

     ps2

 ps1,dir_paperfigs+'mass_colour_data.ps'
  !p.multi=0
  plotsym,0,0.2,/fill
  plot, lgm[ind].per50, uvj[0,ind]-uvj[1,ind],psym=3,ytitle='U-V',xtitle=textoidl('log M*/M_\odot'),xr=[8.7,11.5],yr=[0,2.5],/nodata,/ys,/xs

  vwplot_2dcontour, lgm[ind_red].per50,uvj[0,ind_red]-uvj[1,ind_red],[8.5,12,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[0]),c_linestyle=0
;  oplot, lgm[ind_red].per50,uvj[0,ind_red]-uvj[1,ind_red],psym=8,color=cgcolor(stack_color[0])

  vwplot_2dcontour, lgm[ind_sf2].per50,uvj[0,ind_sf2]-uvj[1,ind_sf2],[8.5,12,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[2]),c_linestyle=2
; oplot,lgm[ind_sf2].per50, uvj[0,ind_sf2]-uvj[1,ind_sf2],psym=8,color=cgcolor(stack_color[2])
 
  vwplot_2dcontour, lgm[ind_sb].per50,uvj[0,ind_sb]-uvj[1,ind_sb],[8.5,12,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[3]),c_linestyle=3
;  oplot, lgm[ind_sb].per50,uvj[0,ind_sb]-uvj[1,ind_sb],psym=8,color=cgcolor(stack_color[3])

  vwplot_2dcontour, lgm[ind_sf3].per50,uvj[0,ind_sf3]-uvj[1,ind_sf3],[8.5,12,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[1]),c_linestyle=4
;  oplot, lgm[ind_sf3].per50,uvj[0,ind_sf3]-uvj[1,ind_sf3],psym=8,color=cgcolor(stack_color[1])

  vwplot_2dcontour, lgm[ind_dusty].per50,uvj[0,ind_dusty]-uvj[1,ind_dusty],[8.5,12,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[5]),c_linestyle=5
;  oplot, lgm[ind_dusty].per50,uvj[0,ind_dusty]-uvj[1,ind_dusty],psym=8,color=cgcolor(stack_color[5])

  plotsym,0,0.4 & oplot,lgm[ind_psb].per50, uvj[0,ind_psb]-uvj[1,ind_psb],psym=8,color=cgcolor('black')
  plotsym,8,0.4 & oplot,lgm[ind_lomet].per50, uvj[0,ind_lomet]-uvj[1,ind_lomet],psym=8,color=cgcolor('black')
  plotsym,0,0.3,/fill & oplot,lgm[ind_psb].per50, uvj[0,ind_psb]-uvj[1,ind_psb],psym=8,color=cgcolor(stack_color[4])
  plotsym,8,0.3,/fill & oplot,lgm[ind_lomet].per50, uvj[0,ind_lomet]-uvj[1,ind_lomet],psym=8,color=cgcolor(stack_color[6])

  x = findgen(10)-23
  y = -0.1*x-0.7
  oplot, x,y,linestyle=2

  y =-0.13*x-1.7
  oplot, x,y,linestyle=2

  ps2


;;------------------------------------------------------------------
;;-- UVJ colours: for comparison with other work
;;------------------------------------------------------------------


  ps1,dir_paperfigs+'uvj_data.ps'
  !p.multi=0
  plotsym,0,0.2,/fill
  plot, uvj[1,ind]-uvj[2,ind],uvj[0,ind]-uvj[1,ind],psym=3,xtitle='V-J',ytitle='U-V',yr=[0.,2.5],xr=[-0.25,2],/nodata,/xs,/ys
 
  vwplot_2dcontour,  uvj[1,ind_red]-uvj[2,ind_red],uvj[0,ind_red]-uvj[1,ind_red],[-0.5,2,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[0]),c_linestyle=0
;  oplot, uvj[1,ind_red]-uvj[2,ind_red],uvj[0,ind_red]-uvj[1,ind_red],psym=8,color=cgcolor(stack_color[0])

  vwplot_2dcontour,  uvj[1,ind_sf2]-uvj[2,ind_sf2],uvj[0,ind_sf2]-uvj[1,ind_sf2],[-0.5,2,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[2]),c_linestyle=2
;  oplot, uvj[1,ind_sf2]-uvj[2,ind_sf2],uvj[0,ind_sf2]-uvj[1,ind_sf2],psym=8,color=cgcolor(stack_color[2])

  vwplot_2dcontour,  uvj[1,ind_sb]-uvj[2,ind_sb],uvj[0,ind_sb]-uvj[1,ind_sb],[-0.5,2,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[3]),c_linestyle=3
;  oplot, uvj[1,ind_sb]-uvj[2,ind_sb],uvj[0,ind_sb]-uvj[1,ind_sb],psym=8,color=cgcolor(stack_color[3])
 
  vwplot_2dcontour,  uvj[1,ind_sf3]-uvj[2,ind_sf3],uvj[0,ind_sf3]-uvj[1,ind_sf3],[-0.5,2,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[1]),c_linestyle=4
;  oplot, uvj[1,ind_sf3]-uvj[2,ind_sf3],uvj[0,ind_sf3]-uvj[1,ind_sf3],psym=8,color=cgcolor(stack_color[1])

  vwplot_2dcontour,  uvj[1,ind_dusty]-uvj[2,ind_dusty],uvj[0,ind_dusty]-uvj[1,ind_dusty],[-0.5,2,0,2.5],[0.1,0.1],psym=8,/log,/overplot,levels=[0.5,0.75,1.,1.25,1.5,1.75,2.0],color=cgcolor(stack_color[5]),c_linestyle=5
;  oplot, uvj[1,ind_dusty]-uvj[2,ind_dusty],uvj[0,ind_dusty]-uvj[1,ind_dusty],psym=8,color=cgcolor(stack_color[5])

  plotsym,0,0.4
  oplot, uvj[1,ind_psb]-uvj[2,ind_psb],uvj[0,ind_psb]-uvj[1,ind_psb],psym=8,color=cgcolor('black')
  plotsym,8,0.4
  oplot, uvj[1,ind_lomet]-uvj[2,ind_lomet],uvj[0,ind_lomet]-uvj[1,ind_lomet],psym=8,color=cgcolor('black')

  plotsym,0,0.3,/fill
  oplot, uvj[1,ind_psb]-uvj[2,ind_psb],uvj[0,ind_psb]-uvj[1,ind_psb],psym=8,color=cgcolor(stack_color[4])
  plotsym,8,0.3,/fill
  oplot, uvj[1,ind_lomet]-uvj[2,ind_lomet],uvj[0,ind_lomet]-uvj[1,ind_lomet],psym=8,color=cgcolor(stack_color[6])


  oplot, [-0.2,-0.05], [2.35,2.35],color=cgcolor(stack_color[3]),linestyle=2 & xyouts, 0, 2.33 ,'SF1',color=cgcolor(stack_color[3])
  oplot, [-0.2,-0.05], [2.25,2.25],color=cgcolor(stack_color[2]),linestyle=4 & xyouts, 0, 2.23 ,'SF2',color=cgcolor(stack_color[2])
  oplot, [-0.2,-0.05], [2.15,2.15],color=cgcolor(stack_color[1]),linestyle=3 & xyouts, 0, 2.13 ,'SF3',color=cgcolor(stack_color[1])
  oplot, [-0.2,-0.05], [2.05,2.05],color=cgcolor(stack_color[0]),linestyle=0 & xyouts, 0, 2.03 ,'Red',color=cgcolor(stack_color[0])
  oplot, [-0.2,-0.05], [1.95,1.95],color=cgcolor(stack_color[5]),linestyle=5 & xyouts, 0, 1.93 ,'Dusty',color=cgcolor(stack_color[5])
  plotsym,0,0.4 & plots, -0.1, 1.85,psym=8,color=cgcolor('black') 
  plotsym,0,0.3,/fill & plots, -0.1, 1.85,psym=8,color=cgcolor(stack_color[4]) 
  xyouts, 0, 1.83,'PSB',color=cgcolor(stack_color[4])
  plotsym,8,0.4 & plots, -0.1,1.75,psym=8,color=cgcolor('black')
  plotsym,8,0.3,/fill & plots, -0.1, 1.75,psym=8,color=cgcolor(stack_color[6]) 
  xyouts, 0., 1.73,'Low-met',color=cgcolor(stack_color[6])


  x = 0.86*findgen(10)/10.+0.74
  y = 0.8*x+0.7
  oplot, x,y 
  oplot, [1.5,1.5],[1.9,2.5]
  oplot, [-.5,0.74],[1.3,1.3]
  oplot, [0.95,0.95],[1.45,2.5],linestyle=1
  ps2


  endif




END
