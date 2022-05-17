;; f_nu in Jansky
;; mag is AB

;; gives same answer for 48.574
;; magAB = -2.5 x log(Fν) + 8.926
;; Fv = 10d^(-0.4*(magAB-8.926))

FUNCTION MAGTOJANSKY, mag_AB

f_nu_ergs = 10d^(-0.4*(mag_AB + 48.6))

;; 1Jy = 10-23 erg s-1 Hz-1 cm-2
f_nu = f_nu_ergs*1d23

return, f_nu


END
