
;  stack_name = ['Red','SF3','SF2','SF1','PSB', 'Dusty','low-Z']  

;  ind = where(strtrim(data.CLASS) eq 'Red') & if ind[0] ne -1 then struc_topcat[ind].CLASSnb = 1
;  ind = where(strtrim(data.CLASS) eq 'SF3') & if ind[0] ne -1 then struc_topcat[ind].CLASSnb = 2
;  ind = where(strtrim(data.CLASS) eq 'SF2') & if ind[0] ne -1 then struc_topcat[ind].CLASSnb = 3
;  ind = where(strtrim(data.CLASS) eq 'SF1') & if ind[0] ne -1 then  struc_topcat[ind].CLASSnb = 4
;  ind = where(strtrim(data.CLASS) eq 'PSB') & if ind[0] ne -1 then struc_topcat[ind].CLASSnb = 5
;  ind = where(strtrim(data.CLASS) eq 'Dusty') & if ind[0] ne -1 then struc_topcat[ind].CLASSnb = 6
;  ind = where(strtrim(data.CLASS) eq 'low-Z' or strtrim(data.CLASS) eq 'Low-met') & if ind[0] ne -1 then struc_topcat[ind].CLASSnb = 7
;  ind = where(strtrim(data.CLASS) eq 'NOCLASS') & if ind[0] ne -1 then struc_topcat[ind].CLASSnb = -1


PRO VWSC_OUTPUT, pcs_data, pcerr_data, pcs_phot, pcerr_phot, ngood, norm, chisq, lines,  fluxfile, massfile, magfile, outfile

  
ngal = n_elements(norm)
namplitudes = (size(pcs_data,/dim))[0]
;;-- Classes, for easy use
VWSC_CLASSIFY, pcs_data, lines,ind_red,ind_sf2, ind_sf1, ind_psb, ind_sf3, ind_dusty, ind_lomet,ind_junk, class, stack_name, stack_color

class = fltarr(ngal)
if ind_red[0] ne -1 then class[ind_red] = 1
if ind_sf3[0] ne -1 then class[ind_sf3] = 2
if ind_sf2[0] ne -1 then class[ind_sf2] = 2
if ind_sf1[0] ne -1 then class[ind_sf1] = 2
if ind_psb[0] ne -1 then class[ind_psb] = 5
if ind_dusty[0] ne -1 then class[ind_dusty] = 6
if ind_lomet[0] ne -1 then class[ind_lomet] = 7
if ind_junk[0] ne -1 then class[ind_junk] = -1


;;-- Catalogue contains flux, flux_err, zphot, ID, [mass, ra, dec], kmag
restore, fluxfile
nband = (size(flux,/dim))[0]

;;-- Contain other parameters
if file_test(massfile) ne 0 then restore, massfile
if file_test(magfile) ne 0 then restore, magfile
if file_test(magfile) ne 0 then uvj = transpose([[ABmag_u_z0],[ABmag_V_z0],[ABmag_J_z0]])
  
;;-- output strucutre for fits format

output = replicate({flux:fltarr(nband), flux_err:fltarr(nband), zphot:0.0, zspec:0.0, zused: 0.0, SC:fltarr(namplitudes), SCerr:fltarr(namplitudes), SCphot:fltarr(namplitudes), SCerrphot:fltarr(namplitudes), norm:0.0, chisq:0.0, ngood:0, ID:0L, CLASS:0L, Kmag:0.0,mass_SC:fltarr(3),lgage_wr_SC:fltarr(3),f_burst_SC:fltarr(3),tau_av_SC:fltarr(3),zmet_SC:fltarr(3),lgsfr_SC:fltarr(3),lgssfr_SC:fltarr(3), ra:0.0, dec:0.0, uvj:fltarr(3)},ngal)

;;changed CLASS from string to integer
if n_elements(Kmag) ne 0 then output.Kmag = Kmag
output.flux = flux
output.flux_err = flux_err
output.zphot = zphot
if n_elements(zspec) ne 0 then output.zspec = zspec
if n_elements(zused) ne 0 then output.zused = zused
output.SC = pcs_data
output.SCerr = pcerr_data
if n_elements(pcs_phot) ne 0 then output.SCphot = pcs_phot
if n_elements(pcs_phot) ne 0 then output.SCerrphot = pcerr_phot
output.norm = norm
output.chisq = chisq
output.ngood = ngood
output.ID = ID
output.CLASS = CLASS


if n_elements(lgm) ne 0 then begin
   output.mass_SC[0] = lgm.per16
   output.mass_SC[1] = lgm.per50
   output.mass_SC[2] = lgm.per84
endif

if n_elements(lgssfr) ne 0 then begin
   output.lgssfr_SC[0] = lgssfr.per16
   output.lgssfr_SC[1] = lgssfr.per50
   output.lgssfr_SC[2] = lgssfr.per84
endif
if n_elements(lgsfr) ne 0 then  begin
   output.lgsfr_SC[0] = lgsfr.per16
   output.lgsfr_SC[1] = lgsfr.per50
   output.lgsfr_SC[2] = lgsfr.per84
endif
if n_elements(lgage_wr) ne 0 then  begin
   output.lgage_wr_SC[0] = lgage_wr.per16
   output.lgage_wr_SC[1] = lgage_wr.per50
   output.lgage_wr_SC[2] = lgage_wr.per84
endif
if n_elements(f_burst) ne 0 then  begin
   output.f_burst_SC[0] = f_burst.per16
   output.f_burst_SC[1] = f_burst.per50
   output.f_burst_SC[2] = f_burst.per84
endif
if n_elements(tau_av) ne 0 then begin
    output.tau_av_SC[0] = tau_av.per16
   output.tau_av_SC[1] = tau_av.per50
   output.tau_av_SC[2] = tau_av.per84
endif
if n_elements(z_met) ne 0 then  begin
   output.zmet_SC[0] = z_met.per16
   output.zmet_SC[1] = z_met.per50
   output.zmet_SC[2] = z_met.per84
endif

if n_elements(ra) ne 0 then output.ra = ra
if n_elements(dec) ne 0 then output.dec = dec
if n_elements(uvj) ne 0 then output.uvj = uvj

;;-- save

mwrfits, output, outfile, /create

END
