pro sfbinfit_go, imglist, noiselist, $        ; inputs
                 fdir = fdir, $
                 xyfile = xyfile, goodflag = goodflag, $
                 dxy_sec = dxy_sec, $
                 zoom = zoom,$
                 nonoisefile = nonoisefile, $
                 background_const = background_const, $ ; use constant background
                 subtract_bkgd = subtract_bkgd, $
                 psfsize = psfsize, $
                 outstr = outstr0, $
                 niter_max = NITER_MAX, $
                 postol=postol, flxtol=flxtol, $
                 qzap = qzap, noqzapcheck = noqzapcheck, $
                 nocheck = nocheck, $
                 nozoom = nozoom

;+
; simple wrapper for running SFBINFIT.PRO on multiple images.
; derived from 08/03/07 version of BINFIT_GO.PRO
; 
; 'imglist' can either be a 3-d array of images or a file with a list of FITS files
; 
; Written by M. Liu (IfA/Hawaii) 08/03/07
;-

if not(keyword_set(ZOOM)) then ZOOM = 100
if not(keyword_set(fdir)) then fdir = ''

BADVAL = -1e6
if n_params() lt 1 then begin
   print, 'pro sfbinfit_go, imglist, (noiselist), [fdir=],'
   print, '               [xyfile=], [goodflag],'
   print, '               [dxy_sec=],'
   print, '               [zoom='+strc(ZOOM)+'],'
   print, '               [outstr=], [nonoisefile], [psfsize=], [niter=], [nocheck]' 
   return
endif


; load images and set up output prefixes
loadct, 15, /silent
print
if (size(imglist, /tname) eq 'STRING') then begin
   loadimf, imglist, imgs, /coadd, fdir = fdir, /silent
   nimgs = (size(imgs))(3)
   readcol2, imglist, ilist, form = 'a', /silent
   infile = imglist  ; name of file with input list
   filebreak, ilist, nvfile = outstr
   if keyword_set(outstr0) then $
      outstr = outstr0+'.'+outstr
   message, 'getting images from input file "'+imglist+'"', /info
endif else begin
   imgs = imglist
   nimgs = n_elements(imgs[0,0,*]);(size(imgs))(3)
   ilist = strc(indgen(nimgs))
   infile = 'passed from array'
   message, 'using images passed from 3-d array', /info
   outstr = outstr0+'_'+ilist
endelse
message, '  number of images to process = '+strc(nimgs), /info
print

; load noise images, if passed
; note: /nonoisefile flag takes precedent, even if noise images are passed
if (n_params() eq 2) and not(keyword_set(nonoisefile)) then begin
   if (size(noiselist, /tname) eq 'STRING') then begin
      readcol2, noiselist, nlist, form = 'a'
      noisefiles = nlist
   endif else $
      noisefiles = noiselist
   nonoise = 0
endif else begin
   noisefiles = strarr(nimgs)
   nonoise = 1
endelse


; load (x,y) position of primary, if specified
; if /goodflag, it means use the info about which images are good in the shifts file
;   assumed to be kept in the 4th column of the shifts file
; all the images are processed by SFBINFIT (even bad ones), 
;   but only the good ones are used to compute the average & std dev   
if keyword_set(xyfile) then begin
   if not(keyword_set(goodflag)) then begin
      readcol2, xyfile, dx, dy, /silent 
      wgood = intarr(nimgs)+1
   endif else $
      readcol2, xyfile, dx, dy, ii, wgood, /silent 
   shiftfile = xyfile
endif else begin
   dx = fltarr(nimgs)
   dy = fltarr(nimgs)
   dxy_sec = [0, 0]
   wgood = intarr(nimgs)+1
   shiftfile = ''
endelse
wgood = round(wgood)


; initialize output variables
sep = fltarr(nimgs)  &  pa = fltarr(nimgs)  &  fratio = fltarr(nimgs)  &  niter = fltarr(nimgs)
x_0 = fltarr(nimgs)  &  x_1 = fltarr(nimgs)  &  y_0 = fltarr(nimgs)  &  y_1 = fltarr(nimgs)
f_0 = fltarr(nimgs)  &  f_1 = fltarr(nimgs)

; loop over images
;wclear
for i=0, nimgs-1 do begin
    print
    print, '****************************************'
    print, 'image ', strc(i+1), ' out of ', strc(nimgs)
    print, '  ', ilist(i)
    print, '****************************************'

    if wgood[i] eq 1 then begin
       
       ; get images
       im0 = imgs(*, *, i)

       ; fix bad pixels
       fixpix, im0, im0 ne BADVAL, imfix, /silent

       ; if shift file passed, compute coordinates for primary and
       ; secondary
       if keyword_set(xyfile) then begin
          xy1 = [dx(i), dy(i)]
          if n_elements(dxy_sec) eq 2 then $
             xy2 = [dx(i)+dxy_sec(0), dy(i)+dxy_sec(1)] $
          else $
             xy2 = [dx(i)+dxy_sec(i,0), dy(i)+dxy_sec(i,1)]
          message, 'using user passed positions', /info
          print, '    xy1 = ', commalist(xy1)
          print, '    xy2 = ', commalist(xy2)
       endif

       ; compute constant background
       if keyword_set(background_const) then begin
          resistant_mean, imfix, 2.5, background_const
          print,' using constant background: '+strc(background_const)
       endif
       
                                ; do fitting
       if keyword_set(xyfile) then begin
          sfbinfit, imfix, $
                    fdir+noisefiles(i), nonoisefile = nonoise, $
                    ss, pp, rr, nn, xx, yy, ff, $          
                    psf, syn, $ ; results from final iteration
                    xy1 = xy1, xy2 = xy2, $
                    zoom = ZOOM, $
                    background_const=background_const, $
                    subtract_bkgd = subtract_bkgd, $
                    psfsize = psfsize, $
                    outstr = outstr(i), $
                    nsave_max = NSAVE_MAX, $
                    niter_max = NITER_MAX, $
                    postol=postol, flxtol=flxtol, $
                    nocheck = nocheck, $
                    lessdisplay = lessdisplay, $
                    qzap = qzap, noqzapcheck = noqzapcheck, $
                    nosave = nosave, $
                    nozoom = nozoom 
       endif else $
          sfbinfit, imfix, $
                    fdir+noisefiles(i), nonoisefile = nonoise, $
                    ss, pp, rr, nn, xx, yy, ff, $          
                    psf, syn, $ ; results from final iteration
                    zoom = ZOOM, $
                    background_const=background_const, $
                    subtract_bkgd = subtract_bkgd, $
                    psfsize = psfsize, $
                    outstr = outstr(i), $
                    niter_max = NITER_MAX, $
                    postol=postol, flxtol=flxtol, $
                    nsave_max = NSAVE_MAX, $
                    nocheck = nocheck, $        
                    lessdisplay = lessdisplay, $
                    qzap = qzap, noqzapcheck = noqzapcheck, $
                    nosave = nosave, $
                    nozoom = nozoom    
       ; store results
       sep[i] = ss  &  pa[i] = pp  &  fratio[i] = rr  &  niter[i] = nn
       x_0[i] = xx[0]  &  x_1[i] = xx[1]  &  y_0[i] = yy[0]  &  y_1[i] = yy[1]
       f_0[i] = ff[0]  &  f_1[i] = ff[1]
       
       ; clear variables, to guard against failed fits
       ss = 0  &  pp = 0  &  rr = 0  &  nn = 0 
       xx = 0  &  yy = 0  &  ff = 0  &  psf = 0  &  syn = 0
       
       if not(keyword_set(nocheck)) then $
          checkcontinue

    endif

endfor


; write out the results
; open output summary file
if not(keyword_set(outstr0)) then $
   outfile = 'sfbinfit_go.RESULTS' $
else $
   outfile = outstr0+'.RESULTS' 
if filecheck(outfile, /over, nocheck = nocheck) then begin
    width = max(strlen(ilist))+200
    openw, unit0, outfile, /get_lun, width = width
    printf, unit0, '#----------------------------------------------------------------------'
    printf, unit0, '# SFBINFIT_GO.PRO: '+systime()+userid()
    printf, unit0, '#   imglist = ', infile
    printf, unit0, '#   shiftfile = ', shiftfile
    printf, unit0, '#   fdir = ', fdir
    if keyword_set(background_const) then $
       printf, unit0, '#   background_const = ', background_const $
    else $
       printf, unit0, '#   background fitted (classic sfbinfit)'
    printf, unit0, '#'
    printf, unit0, '#   ZOOM = ', strc(ZOOM)    
    printf, unit0, '#----------------------------------------------------------------------'
    printf, unit0, '# img   sep(pix)   instrPA    f_0/f_1  niter  good'+$
                   '       x_0           y_0           x_1           y_1          f_0          f_1'
    close, unit0
endif else $
  stop

wg = where(wgood eq 1)
forprint2, out=outfile, /update, width = width, $
           ilist, sep, pa, fratio, niter, round(wgood), x_0, y_0, x_1, y_1, f_0, f_1, $
           f='(a,f10.5,f10.3,f11.5,i6,i6,4f14.6,2e13.3)'
openu, unit0, outfile, /get_lun, /append
printf, unit0, '  '
printf, unit0, 'sep(avg/stddev)     ', string(avg(sep(wg)), '(f9.3)'), $
        '  ', string(stddev(sep(wg)), '(f7.3)')
printf, unit0, 'instPA(avg/stddev)  ', string(avg(pa(wg)), '(f9.3)'), $
        '  ', string(stddev(pa(wg)), '(f7.3)')
printf, unit0, 'f.ratio(avg/stddev) ', string(avg(fratio(wg)), '(f9.3)'), $
        '  ', string(stddev(fratio(wg)), '(f7.3)')
printf, unit0, 'Ntot/Ngood           ', strc(nimgs), '    ', strc(round(total(wgood)))
free_lun, unit0

; show results on screen
spawn, 'more '+outfile
print, '*** FINISHED ***'
end
