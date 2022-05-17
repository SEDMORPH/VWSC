;;*** READ_BC03FILES
;;
;; AIM: read stochastic burst files that have been output from
;; read_spectra*.f and output them into IDL structures and
;; arrays. From here they can be saved in sensible formats as
;; required.  
;;
;; INPUT: infile e.g 'hyb_b' or 'scratch/bc03_lib/rand'
;;
;; OPTIONAL INTPUT: nfiles = nb of files 
;; (the stochastic burst files come in 4 sections, labelled 1-4) if
;; nfiles isn't given then it is assumed that all 4 are required
;; 
;; KEYWORDS: unatten = output unattenuated specta
;; 
;; OUTPUT: params
;;         specarr
;;         wave
;; 
;;------------------------------------------------------------------

PRO vwsc_read_stoch, infile, params, specarr, wave,nfile=nfile,unatten=unatten,maxgal=maxgal,readsfrav=readsfrav


;; don't go much above 20000 with whole array, if you have 4G of RAM!
if n_elements(maxgal) eq 0 then maxgal = 50000                  ;maximum nb of objects in array

;;------------
;; read in data arrays
npix = 6900
if n_elements(nfile) eq 0 then nfile = 4

specarr = dblarr(npix, maxgal)
ngal = 0L
ngal_sep = lonarr(nfile)
for i=0,nfile-1 do begin
    if NOT(KEYWORD_SET(unatten)) then openr,lun,infile+string(i+1,form='(I1.1)')+'.spec',/get_lun $
    else  openr,lun,infile+string(i+1,form='(I1.1)')+'.specunatt',/get_lun
    wave = fltarr(npix)
    readf,lun, wave
    while NOT(EOF(lun)) do begin
        tmp = dblarr(npix)
        readf,lun, tmp
        specarr[*,ngal] = tmp
        ngal = ngal+1
        if ngal eq maxgal then begin
            print, 'VWSC_READ_STOCH WARNING: specarr size too small, returning '+string(maxgal,form='(I0)')+' objects'
            close,lun & free_lun,lun
            if i ne 0 then ngal_sep[i] = ngal-total(ngal_sep[0:i-1]) $
            else ngal_sep[i] = ngal 
            goto, outofloop
        endif
    endwhile
    close, lun & free_lun,lun

    if i ne 0 then ngal_sep[i] = ngal-total(ngal_sep[0:i-1]) $
    else ngal_sep[i] = ngal 

endfor

outofloop:

specarr = specarr[*,0:ngal-1]

if not(keyword_Set(readsfrav)) then params = replicate({modelID:0L, tform:0D, gamma:0D, zmet:0D, $
                    tauv0:0D, mu:0D, nburst:0, mstr1:0D, mstr0:0D, $
                    tlastburst:0D,fburst:dblarr(5), ftot:dblarr(5), $
                    age_wm:0D, age_wr:0D,Dn4000:0D,HdeltaA:0D,HgammaA:0D,Lha:0D,Lhb:0D},ngal)

if (keyword_Set(readsfrav)) then params = replicate({modelID:0L, tform:0D, gamma:0D, zmet:0D, $
                         tauv0:0D, mu:0D, nburst:0, mstr1:0D, mstr0:0D, $
                         tlastburst:0D,fburst:dblarr(5), ftot:dblarr(5), $
                         age_wm:0D, age_wr:0D,Dn4000:0D,HdeltaA:0D,HgammaA:0D,Lha:0D,Lhb:0D,sfrav:dblarr(5)},ngal)

;!!!! if this doesn't work it is because Lha and Lhb have been added
;to some files later (random ones).  sfrav has been added to "all" files, but not
;"rand" files. 
for i=0,nfile-1 do begin
    if ngal_sep[i] eq 0 then continue
    tmp = params[0:ngal_sep[i]-1]
    openr,lun,infile+string(i+1,form='(I1.1)')+'.param',/get_lun
    readf,lun,tmp
    close,lun & free_lun,lun
    if i eq 0 then params[0:ngal_sep[i]-1] = tmp
    if i ne 0 then begin
       params[total(ngal_sep[0:i-1]):total(ngal_sep[0:i-1])+ngal_sep[i]-1] = tmp
       ;; vw bug fixed 24/09 to make sure all the galaxies have
       ;; different modelIDs when using >1 file
       params[total(ngal_sep[0:i-1]):total(ngal_sep[0:i-1])+ngal_sep[i]-1].modelid = params[total(ngal_sep[0:i-1]):total(ngal_sep[0:i-1])+ngal_sep[i]-1].modelid+total(ngal_sep[0:i-1])
    endif
endfor

end
