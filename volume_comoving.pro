;;*** VOLUME_COMOVING
;;
;; AIM: calculate comoving volume between two redshift points
;;      e.g. http://www.astro.ufl.edu/~guzman/ast7939/projects/project01.html 
;;      or Hogg 2000: astro-ph/9905116
;;******************************************************************

FUNCTION comoving, z

common cosm,Olam,Omat,H0 

Otot = Omat+Olam
c = 3E5                     ;km/s

Dl = lumdist(z,H0=H0,Omega_m=Omat, Lambda0=Olam,/silent) ;Mpc

Rodr = (c/H0)*((1-Otot)*(1+z)^2+Olam+Omat*(1+z)^3)^(-0.5)

dV = Rodr*4*!PI*(Dl/(1.+z))^2

return,dV

END

;;------------------------------------------------------------------

FUNCTION volume_comoving, z1, z2,area=area, Olam=Olam_in, Omat=Omat_in,H0=H0_in

common cosm, Olam,Omat,H0

if n_elements(z1) eq 0 or n_elements(z2) eq 0 then begin
    message,'usage: vol=volume_comoving(z1,z2,[area=]) where z1<z2, area=radians'
    return,0
endif

if z1 gt z2 then begin
    message,'usage: vol=volume_comoving(z1,z2,[area=]) where z1<z2, area=radians'
    return,0
endif

if keyword_set(area) then begin
    if area gt !PI*4 then begin
        message,'usage: vol=volume_comoving(z1,z2,[area=]) where z1<z2, area=radians'
        message,'suspicious: area appears to be > 4*PI'
        return,0
    endif
endif

;; default cosmological parameters 
if not(keyword_set(Olam_in)) then Olam=0.7 else Olam=Olam_in
if not(keyword_set(Omat_in)) then Omat=0.3 else Omat=Omat_in
if not(keyword_set(H0_in)) then H0=70. else H0=H0_in


;; return integral of dV = (c/H0) * Dl^2 / (E(z) * (1+z)^2) dA dz
;; between z1 and z2
;; if area not input than use 4!PI, else calculate for given dA

if not(keyword_set(area)) then return, QROMB('comoving',z1,z2) $
else return, area*QROMB('comoving',z1,z2)/(!PI*4)

END
