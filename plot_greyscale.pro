;*** PLOT_GREYSCALE.PRO
;
; AIM: overplot greyscale squares. Designed for use with hist_2d output
;    
; INPUT: Z    = 2D "height" array
;        X, Y = bottom left position of bin
;
; OPTIONAL INPUT: 
;        XDATA = x data array (for overplotting outliers. use rangedata
;                for more control)
;        YDATA = y data array 
;        RANGEDATA = x and y ranges between which outliers should be
;                    overplotted. If not specified, ! values used 
;
;******************************************************************

PRO PLOT_GREYSCALE, Z, X, Y, xdata=xdata, ydata=ydata, log=log, loadct=loadct,range_out = range_out,range_in = range_in,nobar=nobar,rangedata = rangedata,barformat=barformat,bartitle=bartitle,barposition=barposition,CALIFA=CALIFA,ctreverse=ctreverse,range_color=range_color

if N_ELEMENTS(loadct) eq 0 then loadct, 1 else loadct, loadct ;blue,white scale
nx = n_elements(X)
ny = n_elements(Y)

if (size(z,/dim))[0] ne nx or (size(z,/dim))[1] ne ny then begin
    print,'PLOT_GREYSCALE: Z, X, Y arrays incompatible dimensions'
    return
endif

if keyword_set(log) then begin
    Z2 = alog10(Z) 
    if keyword_set(range_in) then range_in2 = alog10(range_in) else range_in2 = minmax(z2[where(z ne 0.)]) 
endif else begin
    Z2 = Z
    if keyword_set(range_in) then range_in2 = range_in else range_in2 = minmax(z2[where(z ne 0.)])  
endelse   

print, 'Range in height',range_in2

dbinx = X[1]-X[0]
dbiny = Y[1]-Y[0]

color = bytscl(Z2, top=255,min=min(range_in2),max = max(range_in2))
if not(keyword_set(CALIFA)) and not(keyword_set(ctreverse)) then color = 255-color               ;white on outskirts


for i=0L,nx -1 do begin 
    for j=0L,ny-1 do begin
        ;; reset colour for Z=0.
        ;if z[i,j] eq 0.0 then continue
       if z2[i,j] lt min(range_in2) or z2[i,j] gt max(range_in2) then continue 
       if i eq nx-1 and j eq ny-1 then $
          polyfill, [X[i], X[i], X[i]+dbinx, X[i]+dbinx], $
          [Y[j], Y[j]+dbiny, Y[j]+dbiny, Y[j]], color=color[i,j] $
        else if i eq nx-1 then $
          polyfill, [X[i], X[i], X[i]+dbinx, X[i]+dbinx], $
          [Y[j], Y[j+1], Y[j+1], Y[j]], color=color[i,j] $
        else if j eq ny-1 then $
          polyfill, [X[i], X[i], X[i+1], X[i+1]], $
          [Y[j], Y[j]+dbiny, Y[j]+dbiny, Y[j]], color=color[i,j] $
          ;; this is the main line...
        else polyfill, [X[i], X[i], X[i+1], X[i+1]], [Y[j], Y[j+1], Y[j+1], Y[j]], color=color[i,j]

      
    endfor
endfor

if n_elements(barposition) eq 0 then begin
   if min(range_in2) ge 0. then position = [0.88, 0.12, 0.93, 0.92] else position = [0.85, 0.12, 0.90, 0.92]
endif else position=barposition

if n_elements(barformat) eq 0 then barformat='(F0.2)'
if not(keyword_set(nobar)) then begin
    cgcolorbar, range = range_in2,/vertical, /right,/reverse,format = barformat,position=position,charsize = !p.charsize-0.2,divisions=8
    if n_elements(bartitle) ne 0 then xyouts, position[0]+0.025, position[3]+0.02, bartitle,/normal,align=0.5
endif

;; plot the extra data points that are in very underdense regions
if n_elements(xdata) ne 0 and n_elements(ydata) ne 0 then begin

    ind = color(reverse(sort(color))) ;high-low : white - dark
    ind2 = where(ind ne 255,count)
    minlevel = (ind[ind2[0.1*count]])[0] ;bottom 10%

    if n_elements(rangedata) eq 0 then rangedata = [!x.crange[0],!x.crange[1],!y.crange[0],!y.crange[1]]

    if n_elements(psym) eq 0 then psym=8

    for i=0L,n_elements(xdata)-1 do begin

        if where(xdata[i] gt rangedata[0] and xdata[i] lt rangedata[1]) eq -1 or $
          where(ydata[i] gt rangedata[2] and ydata[i] lt rangedata[3]) eq -1 then continue

        ;; which colour bin?
        indx = (where(x gt xdata[i]))[0]-1
        indy = (where(y gt ydata[i]))[0]-1        
                
        ;; if outside colour boundary
        if indx le -1 or indy le -1  then plots, xdata[i],ydata[i],psym=psym $
        else $
          ;; if under certain colour boundary level...        
        if color[indx,indy] gt minlevel then plots, xdata[i],ydata[i],psym=psym
        
    endfor
    
endif

range_color = minmax(color)
range_out = minmax(z[where(z ne 0)])

END
