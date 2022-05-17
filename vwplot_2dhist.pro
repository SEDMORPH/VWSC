PRO VWPLOT_2DHIST,xarr,yarr, range,bin, _extra=extra, weight1=weight,weight2=weight2,range_out=range_out,range_in = range_in,loadct=loadct,log=log,nobar=nobar,xrange=xrange,yrange=yrange,barformat=barformat,position=position,nodata=nodata,total=total,bartitle=bartitle,vmax=vmax,barposition=barposition,levels=levels,path_xy=path_xy,path_info=path_info,markcontour=markcontour,x_out=x,y_out=y,z_out=z2,ctreverse=ctreverse,range_color=range_color,noplot=noplot

!x.style=1                      ;can't remember why this is necessary
!y.style=1                      ;but it is!

if n_elements(loadct) eq 0 then loadct = 0

if n_elements(vmax) eq 0 then vmax = fltarr(n_elements(xarr))+1.

nx = abs(range[1]-range[0])/bin[0]
ny = abs(range[3]-range[2])/bin[1]
x = findgen(nx+1)*bin[0] +range[0] ;bottom left coords
y = findgen(ny+1)*bin[1] +range[2]

if n_elements(xrange) eq 0 then xrange = [range[0],range[1]]
if n_elements(yrange) eq 0 then yrange = [range[2],range[3]]


if n_elements(weight) eq 0 then begin
    z = uintarr(nx,ny)
    for i=0, nx-1 do begin
        for j=0,ny-1 do begin
            ind = where(xarr gt x[i] and xarr lt x[i+1] and yarr gt y[j] and yarr lt y[j+1],count)
            if ind[0] ne -1 then z[i,j] = total(1./vmax[ind])
        endfor
    endfor
endif


if n_elements(weight) ne 0. then begin
    z  = fltarr(nx,ny)
    for i=0,nx-1 do begin
        for j=0,ny-1 do begin
            ind = where(xarr gt x[i] and xarr lt x[i+1] and yarr gt y[j] and yarr lt y[j+1]) 
            if ind[0] ne -1 and  n_elements(weight2) eq 0 and not(keyword_set(total)) then $
              z[i,j] = total(weight[ind]/vmax[ind])/total(1./vmax[ind])
            if ind[0] ne -1 and  n_elements(weight2) eq 0 and keyword_set(total) then $
              z[i,j] = total(weight[ind]/vmax[ind])
            if ind[0] ne -1  and  n_elements(weight2) ne 0 then $
              z[i,j] = total(weight[ind]/vmax[ind])/total(weight2[ind]/vmax[ind])
        endfor
    endfor

endif

x=x[0:nx-1]
y=y[0:ny-1]

if n_elements(position) eq 0 then begin
    if keyword_set(log) then z2=alog10(z) else z2=z
    if min(z2,/nan) lt 0. then position=[[0.15,0.12],[0.82,0.92]] $
    else position=[[0.15,0.12],[0.85,0.92]]
endif

if keyword_set(noplot) then return

;xtickform='(F0.1)',ytickform='(F0.1)',

if not(keyword_set(nobar)) then $
  plot, xarr, yarr, xrange=xrange,yrange=yrange,_extra=extra,xtick_get=xtick,ytick_get=ytick,$
  /nodata, position = position else $
  plot, xarr, yarr, _extra=extra, xrange=xrange,yrange=yrange,xtick_get=xtick,ytick_get=ytick,$
  /nodata



plotsym,0,0.1,/fill
if n_elements(weight) eq 0 and not(keyword_set(nodata)) then plot_greyscale, z,x,y,xdata=xarr,ydata=yarr,rangedata=[xrange[0],xrange[1],yrange[0],yrange[1]],log=log,loadct=loadct,nobar=nobar,barformat=barformat,bartitle=bartitle,barposition=barposition,range_in=range_in,ctreverse=ctreverse,range_out=range_out,range_color=range_color $
else plot_greyscale, z,x,y,loadct=loadct,range_out=range_out, range_in = range_in,log=log,nobar=nobar,barformat=barformat,bartitle=bartitle,barposition=barposition,ctreverse=ctreverse,range_color=range_color

xtick_n = replicate(' ',n_elements(xtick+1))
ytick_n = replicate(' ',n_elements(ytick+1))
axis,xaxis=1,xtickint=xtick[1]-xtick[0],xtickname=xtick_n ;top (without labels)
axis,yaxis=1,ytickint=ytick[1]-ytick[0],ytickname=ytick_n ;right (without labels)

;if not(keyword_set(nobar)) then begin
if not(keyword_set(levels)) then levels=[0.5,1.,1.5,2.,2.5,3.0]
nn=n_elements(levels)+1
if keyword_set(log) then z2 = alog10(z) else z2=z
if keyword_set(log) then begin
   contour, z2,x,y,/overplot,c_colors = replicate(cgcolor('white'),nn),c_labels = uintarr(nn),levels = levels,thick=2
   if n_elements(markcontour) ne 0 then if max(z2) gt markcontour then contour, z2,x,y,/overplot,c_colors = cgcolor('green'),c_labels = 0,levels = markcontour,thick=2, path_xy=path_xy,path_info=path_info,/PATH_DATA_COORDS,closed=0
endif
;endif


;;redraw the axes

;axis,xaxis=0,xtickform='(F0.1)';,charsize=1                    ;bottom
;axis,yaxis=0,ytickform='(F0.1)';,charsize=1                    ;left




END
