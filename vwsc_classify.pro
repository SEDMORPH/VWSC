;;*** VWSC_CLASSIFY
;;
;;** AIM: output galaxy SED classifications based on supercolour
;;        amplitudes and a set of lines in supercolour space
;; 
;;** INPUT:
;;        pcs = supercolours [3,ngal]
;;        lines = structure containing slopes and intercepts of lines of disection
;;
;;** OUTPUT: 
;;        indices of galaxies in each class 
;;        class = class number
;;        
;;** EXAMPLE:
;;        lines = {red_sl:0.4, red_int:3.8, psb_sl:0.18,psb_int:3.5,psb2_sl:-0.4,psb2_int:4,dusty_sl:-0.3, dusty_int:-8,lomet_pc3lo:2.5, lomet_pc1hi:10, green_cut:0 }
;;        VWSC_CLASSIFY, pcs_data, lines,ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty, ind_lomet,ind_junk, class, stack_name, stack_color
;; 
;;
;;** NOTE use of vw_setdifference, to remove interesting features in setdifference

PRO VWSC_CLASSIFY, pcs,lines, ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty,ind_lomet,ind_junk,class,stack_name, stack_color, x_psb1=x_psb,y_psb1=y_psb,x_red=x_red,y_red=y_red,x_dusty=x_dusty,y_dusty=y_dusty,y_psb2=y_psb2,x_psb2=x_psb2, xy_psb=xy_psb, xy_red=xy_red, x2_red=x2_red, y2_red=y2_red,x3_red=x3_red,y3_red=y3_red, silent=silent
     
  ngal = (size(pcs,/dim))[1]

  if tag_exist(lines, 'sb_cut') eq 0 then sb_cut = 20 else sb_cut = lines.sb_cut ;this wasn't included in the UDS work, so made back compatible

  x_red = findgen(25)-35
  y_red = lines.red_sl*x_red+lines.red_int

;  x2_red = -6
;  y2_red = [0, 6]

;  x3_red = [-40, -10]
;  y3_red = -3.371          
     
  x_psb = findgen(40)-10
  y_psb = lines.psb_sl*x_psb+lines.psb_int

  x_dusty = findgen(25)-35
  y_dusty = lines.dusty_sl*x_dusty+lines.dusty_int
 
  x_psb2 = findgen(25)-25
  y_psb2 = lines.psb2_sl*x_psb2+lines.psb2_int
  
;; distance from and along PSB line 
  dist_psb = (xy_psb = fltarr(ngal))
  for i=0L,ngal-1 do begin
     dist_psb[i] = vw_pnt_line(pcs[0:1,i],[min(x_psb),min(y_psb)],[max(x_psb),max(y_psb)],pl)
     xy_psb[i] = pl[0]          ; now distance in x direction ;sqrt((pl[0])^2 + (pl[1])^2)
  endfor
     
;; distance from and along red lines 
  dist_red = (xy_red = fltarr(ngal))
  for i=0L,ngal-1 do begin
     dist_red[i] = vw_pnt_line(pcs[0:1,i],[min(x_red),min(y_red)],[max(x_red),max(y_red)],pl)
     xy_red[i] = pl[0]          ;sqrt((pl[0])^2 + (pl[1])^2)
  endfor

;  dist_red2 = (xy_red2 = fltarr(ngal))
;  for i=0L,ngal-1 do begin
;     dist_red2[i] = vw_pnt_line(pcs[0:1,i],[min(x2_red),min(y2_red)],[max(x2_red),max(y2_red)],pl)
;     xy_red2[i] = pl[0]          ;sqrt((pl[0])^2 + (pl[1])^2)
;  endfor

;  dist_red3 = (xy_red3 = fltarr(ngal))
;  for i=0L,ngal-1 do begin
;     dist_red3[i] = vw_pnt_line(pcs[0:1,i],[min(x3_red),min(y3_red)],[max(x3_red),max(y3_red)],pl)
;     xy_red3[i] = pl[0]          ;sqrt((pl[0])^2 + (pl[1])^2)
;  endfor
 
;; distance from and along dusty line 
  dist_dusty = (xy_dusty = fltarr(ngal))
  for i=0L,ngal-1 do begin
     dist_dusty[i] = vw_pnt_line(pcs[0:1,i],[x_dusty[0],y_dusty[0]],[x_dusty[24],y_dusty[24]],pl)
     xy_dusty[i] = pl[0]          ;sqrt((pl[0])^2 + (pl[1])^2)
  endfor
  
;; distance from and along psb2 line 
  dist_psb2 = (xy_psb2 = fltarr(ngal))
  for i=0L,ngal-1 do begin
     dist_psb2[i] = vw_pnt_line(pcs[0:1,i],[x_psb2[0],y_psb2[0]],[x_psb2[24],y_psb2[24]],pl)
     xy_psb2[i] = pl[0]          ;sqrt((pl[0])^2 + (pl[1])^2)
  endfor
  
  ;;-- Red sequence 
;;  ind_red = where(dist_red lt 0 and dist_psb2 lt 0) ;red and blue
  ind_red = where(dist_red lt 0 and dist_dusty gt 0 and dist_psb2 lt 0)
;;  ind_red = where(dist_red lt 0 and pcs[0,*] lt x2_red and pcs[1,*] gt y3_red and dist_psb2 lt 0) ;red and blue  
  ind_dusty = where(pcs[0,*] lt lines.dust_cut and dist_dusty lt 0)
  ind_dusty = vw_setdifference(ind_dusty,ind_red)
  ind_psb = where(dist_psb lt 0 and dist_psb2 gt 0); and pcs[0,*] lt border_psb) 

  ;;-- Blue cloud
  ind_sf3 =where(pcs[0,*] lt lines.green_cut) 
  ind_sf3 = vw_setdifference(ind_sf3,ind_dusty) ;remove dusty objects, if disjoint return green
  ind_sf3 = vw_setdifference(ind_sf3,ind_red)   ;remove red objects, if disjoint return green
  ind_sf3 = vw_setdifference(ind_sf3,ind_psb)   ;remove psb objects, if disjoint return green

  ind_sf1 = where(pcs[0,*] gt sb_cut)
  ind_sf1 = vw_setdifference(ind_sf1,ind_psb) ;remove psb objects, if disjoint return sb
  ind_sf1 = vw_setdifference(ind_sf1,ind_dusty) ;remove dusty objects, if disjoint return sb

  ind_sf2 = where(pcs[0,*] ge lines.green_cut and pcs[0,*] le sb_cut)
  ind_sf2 = vw_setdifference(ind_sf2,ind_dusty) ;remove dusty objects, if disjoint return blue
  ind_sf2 = vw_setdifference(ind_sf2, ind_psb)   ;remove psb objects, if disjoint return blue
  ind_sf2 = vw_setdifference(ind_sf2, ind_red)   ;remove red objects, if disjoint return blue

lometint = lines.lomet_pc3lo
lometsl = lines.lomet_sl

  ;;-- low metallicity
;  ind_lomet = where(pcs[2,*] gt lines.lomet_pc3lo); and pcs[0,*] lt  lines.lomet_pc1hi)
  ind_lomet = where(pcs[2,*] gt lometsl*pcs[0,*]+lometint); and pcs[0,*] lt  lines.lomet_pc1hi)
  ind_lomet = vw_setdifference(ind_lomet,ind_sf3)                                                  ;remove green objects, if disjoint return lomet
  ind_lomet = vw_setdifference(ind_lomet,ind_sf2)                                                   ;remove blue objects, if disjoint return lomet
  ind_lomet = vw_setdifference(ind_lomet,ind_sf1)                                                     ;remove sb objects, if disjoint return lomet
  ind_lomet = vw_setdifference(ind_lomet,ind_dusty)

  ind_red = vw_setdifference(ind_red,ind_lomet)                                                       ;don't allow red-sequence to be low-met
  ind_psb = vw_setdifference(ind_psb,ind_lomet)                                                       ;don't allow PSBs to be low-met
;  ind_dusty = vw_setdifference(ind_dusty,ind_lomet)                                                   ;for good measure
  ind_dusty = vw_setdifference(ind_dusty,ind_psb)                                                   ;for good measure

  ;; finally remove the rubbish
  ind_junk = where((pcs[0,*] eq 0 and pcs[1,*] eq 0) OR (pcs[0,*] lt -50) OR (pcs[1,*] gt 30) OR (pcs[0,*] gt 225) OR (pcs[1,*] lt -40))

  ind_red = vw_setdifference(ind_red,ind_junk) 
  ind_psb = vw_setdifference(ind_psb,ind_junk)
  ind_sf3 = vw_setdifference(ind_sf3,ind_junk)
  ind_sf2 = vw_setdifference(ind_sf2,ind_junk)
  ind_dusty = vw_setdifference(ind_dusty,ind_junk)
  ind_sf1 =  vw_setdifference(ind_sf1,ind_junk)
  ind_lomet = vw_setdifference(ind_lomet,ind_junk)

  if not(keyword_set(silent)) then splog, 'red, blue, green, psb, sb, dusty, lomet, junk: ', size(ind_red,/dim), size(ind_sf2,/dim), size(ind_sf3,/dim), size(ind_psb,/dim), size(ind_sf1,/dim), size(ind_dusty,/dim),size(ind_lomet,/dim), size(ind_junk,/dim)
  if not(keyword_set(silent)) then splog, 'total: red, blue, green, psb, sb, dusty, lomet, junk: ', size(ind_red,/dim)+size(ind_sf2,/dim)+size(ind_sf3,/dim)+size(ind_psb,/dim)+size(ind_sf1,/dim)+ size(ind_dusty,/dim)+size(ind_lomet,/dim)+size(ind_junk,/dim)
  if not(keyword_set(silent)) then splog, 'total: ', ngal
  
;;-- check they are all classified, and none are classified twice!
  if ngal ne  size(ind_red,/dim)+size(ind_sf2,/dim)+size(ind_sf3,/dim)+size(ind_psb,/dim)+size(ind_sf1,/dim)+ size(ind_dusty,/dim)+size(ind_lomet,/dim)+size(ind_junk,/dim) then stop

  stack_color = ['red','cyan','blue','purple','orange','brown','green']
  stack_name = ['Red','SF3','SF2','SF1','PSB', 'Dusty','low-Z']  
  
  class = uintarr(ngal)
  if ind_red[0] ne -1 then class[ind_red] = 1
  if ind_sf3[0] ne -1 then  class[ind_sf3] = 2
  if ind_sf2[0] ne -1 then  class[ind_sf2] = 2
  if ind_sf1[0] ne -1 then  class[ind_sf1] = 2
  if ind_psb[0] ne -1 then  class[ind_psb] = 5
  if ind_dusty[0] ne -1 then  class[ind_dusty] = 6
  if ind_lomet[0] ne -1 then  class[ind_lomet] = 7
  if ind_junk[0] ne -1 then  class[ind_junk] = -1




END
