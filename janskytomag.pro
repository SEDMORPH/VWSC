;; f_nu_rest in Jansky
;; mag is AB

FUNCTION JANSKYTOMAG, f_nu

;; 1Jy = 10-23 erg s-1 Hz-1 cm-2
;f_nu_ergs = f_nu/1e23
;mag_AB = -alog10(f_nu_ergs)*2.5 - 48.6

mag_AB = -2.5*alog10(f_nu/3631.)

return, mag_AB


END
