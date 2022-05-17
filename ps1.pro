pro ps1,filename,mnras=mnras,_extra=extra

set_plot,'ps',_extra=extra

;;-- this needs to =0 to get nice fonts in mice_figs.pro
!p.font = 0
device,file=filename,/color,bits_per_pixel=8,/bold

;; size for MNRAS columns
if keyword_set(MNRAS) then device,/portrait,xoffset=0,yoffset=0,ysize=6.5,xsize=8.5 else $
   device,/portrait,xoffset=0,yoffset=0,ysize=13,xsize=17
!p.thick=(!x.thick=(!y.thick=4))
!p.charthick=4


end
