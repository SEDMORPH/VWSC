PRO VWSC_PLOT_LFMF, psfile, mag_abs, mass, vmax, vsurvey, ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty,ind_lomet,stack_name, stack_color,mass_title=mass_title, mag_title=mag_title

  if not(keyword_set(mass_title)) then mass_title = 'Log(M*/M_\odot)' ;don't know which mass has been entered
  if not(keyword_set(mag_title)) then mag_title = 'Magnitude'

  ps1,psfile
  device,/portrait,xoffset=0,yoffset=0,ysize=12,xsize=30
  !x.charsize=(!y.charsize=2)
  !y.style=1 & !x.style=1
  !p.multi=[0,3,1]
  !x.margin = [15,3]
  !y.margin=[8,2]

;;-- LF 
  plotsym, 0,/fill 
  plot_histogram,mag_abs,nbins=20,min=-24.5,max=-20,/perunitx,xr=[-20.,-25],ytitle=textoidl('log \phi / Mpc^{-3} dex^{-1}'),xtitle=textoidl(mag_title),weights=1/vmax/vsurvey,vsurvey=vsurvey,compllim=0.5,/ylog,yr=[-6,-2],psym=-8
  plot_histogram,mag_abs[ind_red],nbins=20,min=-24.5,max=-20,color=cgcolor(stack_color[0]),/over,/perunitx,weights=1/vmax[ind_red]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-1,/ylog
  plotsym, 0
  plot_histogram,mag_abs[ind_sf2],nbins=16,min=-24,max=-20,color=cgcolor(stack_color[2]),/over,/perunitx,weights=1/vmax[ind_sf2]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-8,/ylog
  plot_histogram,mag_abs[ind_psb],nbins=10,min=-24.5,max=-20,color=cgcolor(stack_color[4]),/over,/perunitx,weights=1/vmax[ind_psb]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-5,/ylog
  plot_histogram,mag_abs[ind_sf1],nbins=13,min=-22.8,max=-20,color=cgcolor(stack_color[3]),/over,/perunitx,weights=1/vmax[ind_sf1]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-4,/ylog

  plot_histogram,mag_abs[ind_sf3],nbins=20,min=-24.5,max=-20,color=cgcolor(stack_color[1]),/over,/perunitx,weights=1/vmax[ind_sf3]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-2,/ylog
  plot_histogram,mag_abs[ind_dusty],nbins=10,min=-24.5,max=-20,color=cgcolor(stack_color[5]),/over,/perunitx,weights=1/vmax[ind_dusty]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-6,/ylog
  plot_histogram,mag_abs[ind_lomet],nbins=10,min=-24.5,max=-20,color=cgcolor(stack_color[6]),/over,/perunitx,weights=1/vmax[ind_lomet]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-7,/ylog

  psym = -[1,2,8,4,5,6,7]
  for i=0,n_elements(stack_name)-1 do xyouts, 0.28,0.935-i*0.035,stack_name[i],color=cgcolor(stack_color[i]),/norm,charsize=0.8
  for i=0,n_elements(stack_name)-1 do plots, [0.24,0.27],[0.945-i*0.035,0.945-i*0.035], psym=psym[i],color=cgcolor(stack_color[i]),/norm
  

  zz = 1.01                     ;this is the redshift where we want the LF (median(z[ind_loz])
  
  alpha = -1.07                 ;these are Michaele's fitted parameters
  Ms0 = -22.26
  zm = 1.78
  km = 0.47
  zphi = 1.7
  kphi = 1.47
  phi0 = 3.5e-3

  Ms = Ms0 - (zz/zm)^km
  phi = phi0*exp(-(zz/zphi)^kphi)
  xx = -25.5+0.3*findgen(20)
  deltaM = xx - Ms
  oplot, xx, alog10(0.4*alog(10)*phi*10d^(-0.4*deltaM*(alpha+1))*exp(-10d^(-0.4*deltaM))),linestyle=1


;;-- MF
  plotsym, 0,/fill 
  plot_histogram,mass,nbins=25,min=9,max=12.5,/perunitx,xr=[9,12.2],ytitle=textoidl('log \phi / Mpc^{-3} dex^{-1}'),xtitle=textoidl(mass_title), weights=1/vmax/vsurvey,/ylog,yr=[-6,-2],vsurvey=vsurvey,compllim=0.5,psym=-8
  plot_histogram,mass[ind_red],nbins=25,min=9,max=12,color=cgcolor(stack_color[0]),/over,/perunitx,weights=1/vmax[ind_red]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-1,/ylog
  plotsym, 0
  plot_histogram,mass[ind_sf2],nbins=25,min=9,max=12,color=cgcolor(stack_color[2]),/over,/perunitx,weights=1/vmax[ind_sf2]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-8,/ylog
  plot_histogram,mass[ind_psb],nbins=15,min=9,max=12,color=cgcolor(stack_color[4]),/over,/perunitx,weights=1/vmax[ind_psb]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-5,/ylog
  plot_histogram,mass[ind_sf1],nbins=25,min=9,max=12,color=cgcolor(stack_color[3]),/over,/perunitx,weights=1/vmax[ind_sf1]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-4,/ylog
  plot_histogram,mass[ind_sf3],nbins=25,min=9,max=12,color=cgcolor(stack_color[1]),/over,/perunitx,weights=1/vmax[ind_sf3]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-2,/ylog
  plot_histogram,mass[ind_dusty],nbins=15,min=10,max=12,color=cgcolor(stack_color[5]),/over,/perunitx,weights=1/vmax[ind_dusty]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-6,/ylog
  plot_histogram,mass[ind_lomet],nbins=25,min=9,max=12,color=cgcolor(stack_color[6]),/over,/perunitx,weights=1/vmax[ind_lomet]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-7,/ylog


;;-- MF red and blue
  
  plotsym,0,/fill 

  plot_histogram,mass,nbins=25,min=9,max=12.5,/perunitx,xr=[9,12.2],ytitle=textoidl('log \phi / Mpc^{-3} dex^{-1}'),xtitle=textoidl(mass_title),weights=1/vmax/vsurvey,/ylog,yr=[-6,-2],vsurvey=vsurvey,compllim=0.5,psym=-8

  tmp = setunion(setunion(ind_red,ind_psb),ind_lomet)
  plot_histogram,mass[tmp],nbins=25,min=9,max=12.5,color=cgcolor('tomato'),/over,/perunitx,weights=1/vmax[tmp]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-8,/ylog

  tmp = setunion(setunion(setunion(ind_sf1,ind_sf2),ind_sf3),ind_dusty)
  plot_histogram,mass[tmp],nbins=25,min=9,max=12.5,color=cgcolor('navy'),/over,/perunitx,weights=1/vmax[tmp]/vsurvey,vsurvey=vsurvey,compllim=0.5,psym=-8,/ylog

  xyouts, 0.78,0.26-0*0.035,'Total',/norm,charsize=0.8
  plots, [0.74,0.77],[0.27-0*0.035,0.27-0*0.035], psym=-8,/norm
  xyouts, 0.78,0.26-1*0.035,'Red+PSB+low-Z',color=cgcolor('tomato'),/norm,charsize=0.8
  plots, [0.74,0.77],[0.27-1*0.035,0.27-1*0.035], psym=-8,color=cgcolor('tomato'),/norm
  xyouts, 0.78,0.26-2*0.035,'SF1+SF2+SF3+dusty',color=cgcolor('navy'),/norm,charsize=0.8
  plots, [0.74,0.77],[0.27-2*0.035,0.27-2*0.035], psym=-8,color=cgcolor('navy'),/norm
  
;;-- Moustakas et al. 2013 h70, Chabrier IMF , bursty SFH, ELAIS-S1 CDFS XMM-LSS ,
;;   0.8<z<1.0  spectroscopic redshifts, 120,000 galaxies 
  oplot, [10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12.0], [-2.583,-2.658,-2.701,-2.842,-3.039,-3.296,-3.453,-3.77,-4.32,-4.44,-5.07],psym=1
  oplot, [10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8],[-2.639,-2.697,-2.783,-2.873,-2.944,-2.908,-3.011,-3.113,-3.279,-3.405,-3.63,-3.92,-4.22,-5.05,-4.81,-5.42],psym=1,color=cgcolor('blue')
  oplot, [10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8],[-2.861,-2.913,-2.914,-3.040,-3.283,-3.566,-3.63,-3.95,-4.41,-4.68,-5.33],psym=1,color=cgcolor('red')

;;-- ILbert 2013 : full SFing Q
;;-- Chabrier IMF h=0.7, BC03 models, exponentially declining SFH,
;;   expect 0.1dex offset to higher masses, worse in older objects
  Ms = 10d^[10.87,10.81,10.80]
  phi1 = [2.03, 1.57, 0.5]*1e-3
  alpha1 = [-0.52,-0.11,-0.67]
  phi2 = [0.29,0.0,0.48]*1e-3
  alpha2 = [-1.62,0.0,-1.51]
  xx = 10d^(9.5+2*findgen(25)/20.)
  phi = exp(-xx/Ms[0])*(phi1[0]*(xx/Ms[0])^alpha1[0] + phi2[0]*(xx/Ms[0])^alpha2[0])*alog(10)*xx/Ms[0]
  phi_red = exp(-xx/Ms[1])*(phi1[1]*(xx/Ms[1])^alpha1[1] + phi2[1]*(xx/Ms[1])^alpha2[1])*alog(10)*xx/Ms[1]
  phi_sf = exp(-xx/Ms[2])*(phi1[2]*(xx/Ms[2])^alpha1[2] + phi2[2]*(xx/Ms[2])^alpha2[2])*alog(10)*xx/Ms[2]

  oplot,alog10(xx),alog10(phi)
  oplot,alog10(xx),alog10(phi_red),color=cgcolor('red')
  oplot,alog10(xx),alog10(phi_sf),color=cgcolor('blue')

  ps2


END
