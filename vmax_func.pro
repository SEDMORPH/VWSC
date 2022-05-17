;;-- function for vmax_bc03

FUNCTION vmax_FUNC, z

;;-- absmag is the absmag of the galaxy
common func, absmag, appmag, kcorr_igal, zgrid


ind = where(z lt 0.0001)
zin = z
if ind[0] ne -1 then zin[ind]=0.0001

;; seem to have fixed this by increasing tolerance in fx_root
;if size(zin,/type) eq 9 then zin = double(zin)

;lumdist is in Mpc
DM = 5.*alog10(lumdist(zin,/silent)*100000) ;10pc


Kcorr = fltarr(n_elements(z))
for i=0,n_elements(z)-1 do begin
   tmp = zin[i]
   diff = abs(zgrid-tmp)
   ind_grid = where(diff eq min(diff))
   Kcorr[i] = Kcorr_igal[ind_grid]
endfor

;print, ' '
;print, z
;print, kcorr

result = (absmag+DM-Kcorr)-appmag
;print, result

;help, z

return, result


END
