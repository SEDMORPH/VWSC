PRO vwsc_plot_sc123, psfile, modelfile, pcs_data, ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty,ind_lomet,stack_name, stack_color, range12, range13, dbin12, dbin13

  restore, modelfile

  if n_elements(range12) eq 0 then range12 = [-40,70,-15,10]
  if n_elements(range13) eq 0 then range13 = [-40,70,-5,8]
  if n_elements(dbin12) eq 0 then dbin12 = [2.5,0.5]
  if n_elements(dbin13) eq 0 then dbin13 = [2.5,0.25]

  ps1, psfile
  device,/portrait,xoffset=0,yoffset=0,ysize=20,xsize=30
  multiplot, [2,2],ygap=0.01,xgap=0.01,mytitoffset=-1,/doyaxis
  
  ;;-- PC12 histogram
  vwplot_2dhist, pcs_data[0,*],pcs_data[1,*],range12,dbin12,psym=8,xr=range12[0:1],yr=range12[2:3],/xs,/ys,ytitle='SC2',/log,/nobar,levels=[0.75,1.2,1.4,1.6,1.75],loadct=10
  

  ;;-- PC12 models + spectra points 
  multiplot
  vwplot_2dhist, pcs[0,*],pcs[1,*],range12,dbin12,psym=8,xr=range12[0:1],yr=range12[2:3],/xs,/ys,/log,/nobar,loadct=0


 ;;-- PC12 points 
;  multiplot
  plotsym,0,0.1,/fill
;  plot,  pcs_data[0,*],pcs_data[1,*],psym=8,xr=range12[0:1],yr=range12[2:3],/xs,/ys,/nodata
  if ind_red[0] ne -1 then oplot,  pcs_data[0,ind_red],pcs_data[1,ind_red],color=cgcolor(stack_color[0]),psym=8
  if ind_sf2[0] ne -1 then oplot,  pcs_data[0,ind_sf2],pcs_data[1,ind_sf2],color=cgcolor(stack_color[2]),psym=8
  if ind_sf3[0] ne -1 then oplot,  pcs_data[0,ind_sf3],pcs_data[1,ind_sf3],color=cgcolor(stack_color[1]),psym=8
  if ind_sf1[0] ne -1 then oplot,  pcs_data[0,ind_sf1],pcs_data[1,ind_sf1],color=cgcolor(stack_color[3]),psym=8
  if ind_dusty[0] ne -1 then oplot,  pcs_data[0,ind_dusty],pcs_data[1,ind_dusty],color=cgcolor(stack_color[5]),psym=8
;  plotsym,0,0.3,/fill
  if ind_psb[0] ne -1 then oplot,  pcs_data[0,ind_psb],pcs_data[1,ind_psb],color=cgcolor(stack_color[4]),psym=8
  if ind_lomet[0] ne -1 then oplot,  pcs_data[0,ind_lomet],pcs_data[1,ind_lomet],color=cgcolor(stack_color[6]),psym=8


  ;;-- PC13 histogram
  multiplot
  vwplot_2dhist, pcs_data[0,*],pcs_data[2,*],range13,dbin13,psym=8,xr=range13[0:1],yr=range13[2:3],/xs,/ys,ytitle='SC3',xtitle='SC1',/log,/nobar,levels=[0.75,1.25,1.5,1.75],loadct=10

 
  ;;-- PC13 models + spectra points
  multiplot
  vwplot_2dhist, pcs[0,*],pcs[2,*],range13,dbin13,psym=8,xr=range13[0:1],yr=range13[2:3],/xs,/ys,xtitle='SC1',/log,/nobar,loadct=0

 ;;-- PC13 points + spectra points
;  multiplot
  plotsym,0,0.1,/fill
;  plot,  pcs_data[0,*],pcs_data[2,*],psym=3,xr=range13[0:1],yr=range13[2:3],/xs,/ys,xtitle='SC1',/nodata
  if ind_red[0] ne -1 then oplot,  pcs_data[0,ind_red],pcs_data[2,ind_red],color=cgcolor(stack_color[0]),psym=8
  if ind_sf2[0] ne -1 then oplot,  pcs_data[0,ind_sf2],pcs_data[2,ind_sf2],color=cgcolor(stack_color[2]),psym=8
  if ind_sf3[0] ne -1 then oplot,  pcs_data[0,ind_sf3],pcs_data[2,ind_sf3],color=cgcolor(stack_color[1]),psym=8
  if ind_sf1[0] ne -1 then oplot,  pcs_data[0,ind_sf1],pcs_data[2,ind_sf1],color=cgcolor(stack_color[3]),psym=8
  if ind_dusty[0] ne -1 then oplot,  pcs_data[0,ind_dusty],pcs_data[2,ind_dusty],color=cgcolor(stack_color[5]),psym=8
;  plotsym,0,0.3,/fill
  if ind_psb[0] ne -1 then oplot,  pcs_data[0,ind_psb],pcs_data[2,ind_psb],color=cgcolor(stack_color[4]),psym=8
  if ind_lomet[0] ne -1 then oplot,  pcs_data[0,ind_lomet],pcs_data[2,ind_lomet],color=cgcolor(stack_color[6]),psym=8


  multiplot,/reset
  ps2


END
