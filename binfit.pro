pro binfit, image, $
            mask=mask, $ ; 
            sep, pa, fittedparms, model, dmag, $   ; fitting results
            xy1=xy1, xy2=xy2, $              ; initial guesses for (x,y) values
            xyfixed = xyfixed, $             ; use fixed values for (x,y) of primary, don't fit for these
            cbox=CBOX, zoom=ZOOM, fwhm=FWHM0, nozoom = NOZOOM, nofix = nofix, $
            twogauss=TWOGAUSS, threegauss=threegauss, elliptical=elliptical, $
            skyfit=SKYFIT, skyslope=SKYSLOPE, $
            qzap = qzap, nocheckqzap = nocheckqzap, $
            tvscl=tv, $
            noswap = noswap, $
            nodisplay = nodisplay, silent=silent

;+
; given an image of a binary, fit the components with an analytic PSF
; model.  default analytic model is a 2-d circular gaussian + zero sky
; initial guesses for (x,y) positions are made via user interaction.
;
; works well for binaries with modest flux ratios + separations.
; 
; *NOTE: for reasons I cannot figure out, it can be difficult to run
;        this program the first time, b/c IDL doesn't recognize the
;        MPEVALEXPR() function, which is contained within the MPFITEXPR.PRO
;        routine. I seem to be able to get it to work if I manually go ".run
;        MPEVALEXPR", but this has to be done on the copy I have in
;        ~/idl/idp3 and not the copy in ~/idl/markwardt.  When it crashes the
;        first time, do: 
;            IDL> .run mpfitexpr 
;            IDL> .run binfit
; 
; INPUTS
;   image   either IDL 2-d array or name of a FITS file
;   xy1     initial guess for position of primary
;   xy2     initial guess for position of secondary 
; 
; KEYWORD PARMS
;   cbox      centering box size for picking components (default = 7)
;   zoom      size of zoomed region used to pick the secondary (default = 40)
;   fwhm0     initial guess for FWHM (default = 4)
;   twogauss  use  2 gaussian PSF model (narrow & broad component)
;
; USES
;   ireg, gaussian2d (my version), mpfitexpr, binparms, gausspsfint
;
; Written by M. Liu (IfA/UH) 11/16/01 
; 05/06/05 MCL: Makes both linear & log plots of data vs fitted model
; 05/08/05 MCL: Added 'twogauss' capability and 'skyfit' option
;               Final model computation now done with MPEVALEXPRPRO
;               Don't open windows if already open (useful for multiple calls)
;               Very minor bug fix of normlization factor for initial guess of peaks
; 06/13/05 MCL: Added /tv
; 09/03/05 MCL: Added /elliptical to do 2-d elliptical gaussian fitting
;                 but only for case where gaussian is aligned along x&y axes
;                 uses GAUSSPSFINT.PRO
;               Added /silent
; 07/30/06 MCL: Centers the fitted region on the binary midpoint, instead of primary
;               Added plot of minor axis cut of data vs fitted model
;               Added 'threegauss' capability
;               Bug fix for /skyfit (wasn't generating array of initial guesses right)
;               Improved /elliptical to allow for arbitrary orientation, not fixed to x/y axes
; 07/30/06 MCL: now uses (my version of) GAUSSIAN2D.PRO instead of GAUSSPSF.PRO as the core routine
;               now able to request elliptical PSF for all options using /elliptical
;                 previously, /elliptical meant only 2-gaussian elliptical PSF
;                 note: now if set /elliptical only, will mean 1-gaussian elliptical PSF
;               automatically fixes BADVAL pixels
; 08/06/06 MCL: slightly tweaks to displays
; 02/07/07 MCL: changed display range for (quick) primary-subtracted image
; 03/27/07 MCL: added /nozoom (kludge fix for small images)
; 04/27/07 MCL: uses ASINH_STRETCH() for displaying images
; 05/03/07 MCL: added /qzap option
;               added a contour plot that compares data to model 
; 05/15/07 MCL: no longer prints out initial guess of flux ratio (easily confused for final fit)
; 07/02/07 MCL: check if /nozoom should be set (i.e. if image is too small for specified ZOOM value)
;               big fix: was not fixing BADVAL pixels if /nozoom set
;               big fix: was not generating initial guesses if /nozoom set
; 08/05/07 MCL: changed keyword of calls to ROT from cubic=-1 to cubic=-0.5, as recommended by IDL help page
; 08/10/07 MCL: made rotation of model (for plotting) conform to PA convention from BINPARMS.PRO
;               ensured that the primary is forced in first set of variables
;
; --> should allow elliptical PSFs to have different ellipticities
; --> should allow passing of noise map, so can do weighting for fitting
;-

forward_function mpfit
forward_function mpfitexpr

BADVAL = -1e6

if not(keyword_set(CBOX)) then CBOX = 7
if not(keyword_set(ZOOM)) then ZOOM = 100
if not(keyword_set(FWHM0)) then FWHM0 = 4.
if n_params() lt 1 then begin
    print, 'pro binfit, image, sep, pa, fittedparms, '
    print, '            [xy1=], [xy2=], '
    print, '            [cbox=', strc(CBOX), '], [ZOOM = ', strc(ZOOM), '],'+ $
           '            [FWHM0 = '+strc(FWHM0)+'], [twogauss], [threegauss],'
    print, '            [elliptical], [skyfit], [skyslope]'
    return
endif
if keyword_set(twogauss) and keyword_set(threegauss) then begin
    message, '** cannot have /twogauss and /threegauss at same time!**'
    return
endif


; load image
if size(image, /tname) eq 'STRING' then $
 img = readfits(image) $
else $
  img = image
loadct, 15, /silent
sz = size(img)
if ((sz(1) lt ZOOM) or (sz(2) lt ZOOM)) and not(keyword_set(nozoom)) then begin
   message, '** image size is '+strc(sz(1))+' x '+strc(sz(2))+', while ZOOM='+strc(ZOOM)+' **', /info
   message, '   setting /nozoom', /info
   nozoom = 1
endif


; mark initial location, if not already passed
if not(keyword_set(xy1)) then begin
    print, '--------------------------------------------------'
    print, 'STEP 1: mark location of the primary'
    print, '--------------------------------------------------'
    ireg, img, x1, y1, cbox = CBOX, /silent, tv = tv, /asinh
endif else begin
    x1 = xy1(0)
    y1 = xy1(1)
endelse


; extract sub-region of image
if not(keyword_set(nozoom)) then begin
   cut0 = imcut(img-median(img(where(img ne BADVAL))), ZOOM, x1, y1, /silent) 
   x1cut = ZOOM/2 + (x1-round(x1))
   y1cut = ZOOM/2 + (y1-round(y1))
   peak1 = cut0(ZOOM/2, ZOOM/2)
endif else begin
   cut0 = img-median(img(where(img ne BADVAL)))
   x1cut = x1
   y1cut = y1
   peak1 = img(x1cut, y1cut)
endelse
if not keyword_set(nofix) then $
   fixpix, cut0, cut0 ne BADVAL, cut0, /silent


; do crude subtraction of primary via rotation
sub = cut0 - rot(cut0, 180, 1.0, x1cut, y1cut, /pivot, cubic = -0.5)
if not(keyword_set(xy2)) then begin
    print, '----------------------------------------------------'
    print, 'STEP 2: mark the secondary (pick the positive image)'
    print, '----------------------------------------------------'
    ireg, sub, x2cut, y2cut, cbox = CBOX, tv = tv, /silent, /asinh
endif else begin
    x2cut = x1cut + (xy2(0)-xy1(0))
    y2cut = y1cut + (xy2(1)-xy1(1))
 endelse
if round(x2cut) ge n_elements(sub[*,0]) then x2cut=n_elements(sub[*,0])-1
if round(y2cut) ge n_elements(sub[0,*]) then y2cut=n_elements(sub[0,*])-1
peak2 = sub(round(x2cut), round(y2cut))


; extract final sub-region, centered on midpoint of the binary
; remove BADVAL pixels
delta_x = round((x2cut-x1cut)/2.)
delta_y = round((y2cut-y1cut)/2.)
if not(keyword_set(nozoom)) then begin
   ;xoff=x1+delta_x - ZOOM/2
   ;yoff=y1+delta_y - ZOOM/2
   cutnew = imcut(img-median(img(where(img ne BADVAL))), ZOOM, $
                  x1+delta_x, y1+delta_y, xc1, yc1, xoff, yoff, /silent)
   x1 = ZOOM/2 - delta_x
   y1 = ZOOM/2 - delta_y
   x2 = ZOOM/2 + delta_x
   y2 = ZOOM/2 + delta_y
endif else begin
   xoff=0
   yoff=0
   cutnew = img-median(img(where(img ne BADVAL)))
   ;x1 = xy1(0)
   ;y1 = xy1(1)
   ;x2 = xy2(0)
   ;y2 = xy2(1)
   x2 = x2cut;x1 + delta_x
   y2 = y2cut;y1 + delta_y
endelse
if not keyword_set(nofix) then $
   fixpix, cutnew, cutnew ne BADVAL, cut, /silent $
else $
   cut = cutnew


; if desired, run QZAP on the sub-section
; (compute stats on the whole image, since sub-section is usually too small)
if keyword_set(qzap) then begin
   message, 'running QZAPBOTH to remove cosmic rays', /info
   iterstat, img, iout, /silent
   qzapboth, cut, cutfix, sigma = iout(3), /silent;, box=3, nsig=5.
   if not(keyword_set(nocheckqzap)) then begin
      loadct, 15
      win, 0, /checkopen
      display2, asinh_stretch(cut), /tv, tit = 'original'
      win, 1, /checkopen
      display2, asinh_stretch(cutfix), /tv, tit = 'QZAP-ed'
      ;blink, [0, 1]
   endif
   cut = temporary(cutfix)
endif


; if desired, remove median flux level
if keyword_set(medsub) then begin
    message, 'removing median level = '+strc(median(image)), /info
    cut = cut - median(image)
endif


; generate initial guesses
; remember that GAUSSPSF.PRO returns a gaussian which has integrated flux = 1.0
;   thus we express the different components with the integrated flux
;   the peak flux is related by:
;      (integrated) = peak * (2*!pi) * (FWHM/2.35)^2.
ff0 = (2*!pi) * (FWHM0/2.35)^2.    ; conversion from peak flux to integrated flux
;ff0 = sqrt(2*!pi) * FWHM0/2.35    ; original version
if keyword_set(twogauss) then begin
    start = [peak1*ff0, x1, y1, $  ; primary: flux(gaussian-a), x, y
             peak2/peak1, x2, y2, $  ; 2ndary: flux(gaussian-a), x, y
             FWHM0, $              ; FWHM(a)
             2., FWHM0*2.]         ; relative flux of 2nd gaussian, FWHM(b)
             ;0.1, FWHM0*2.]   ; relative *peak* of 2nd gaussian to 1st, FWHM(b)
endif else if keyword_set(threegauss) then begin
    start = [peak1*ff0, x1, y1, $  ; primary: flux(gaussian-a), x, y
             peak2/peak1, x2, y2, $  ; 2ndary: flux(gaussian-a), x, y
             FWHM0, $              ; FWHM(a)
             2., FWHM0*2., $       ; relative flux of 2nd gaussian, FWHM(b)
             4., FWHM0*4.]         ; relative flux of 3rd gaussian, FWHM(c)
             ;0.1, FWHM0*2.]   ; relative *peak* of 2nd gaussian to 1st, FWHM(b)
endif else begin
    start = [peak1*ff0, x1, y1, $    ; primary: flux, x, y
             peak2/peak1, x2, y2, $    ; 2ndary: flux, x, y
             FWHM0]                  ; FWHM
endelse

; allow for elliptical gaussian, if desired
if (keyword_set(elliptical)) then begin
    estr = 'allowing for *elliptical* PSF'
    np = n_elements(start)
    ess = '* p['+strc(n_elements(start))+'], p['+strc(n_elements(start)+1)+']'
    ii_eratio = n_elements(start)
    ii_pa = n_elements(start)+1
    start = [start, $
             1.1, $   ; ratio of FWHM(xaxis) / FWHM(yaxis)
             0.5]     ; initial guess of PA, in radians
endif else begin
    estr = 'assuming circular PSF'
    ess= '* 1.0, 0.0'
endelse


; generate string expression for the analytic model
if not(keyword_set(silent)) then begin
    print, '------------------------------------------------------------'
    print, 'STEP 3: doing fit (be patient)'
    print, '------------------------------------------------------------'
endif
if not(keyword_set(nozoom)) then begin
   ss = strc(ZOOM)+', '+strc(ZOOM)
endif else begin
   sz = size(cut)
   ss = strc(sz(1))+', '+strc(sz(2))
endelse
if (keyword_set(twogauss)) then begin
    ; fit using 2-gaussian PSF, using circular gaussians
    if not(keyword_set(silent)) then $
      print, 'using double-gaussian model for PSF (be extra patient)'
    expr = '     p[0] * gaussian2d('+ss+', p[1], p[2], p[6]/2.35, p[6]/2.35'+ess+', /norm) + '+ $
           'p[3]*p[0] * gaussian2d('+ss+', p[4], p[5], p[6]/2.35, p[6]/2.35'+ess+', /norm) + '+ $
           '     p[0]*p[7] * gaussian2d('+ss+', p[1], p[2], p[8]/2.35, p[8]/2.35'+ess+', /norm) + '+ $
           'p[3]*p[0]*p[7] * gaussian2d('+ss+', p[4], p[5], p[8]/2.35, p[8]/2.35'+ess+', /norm)'
    ;expr = 'p[0] * gausspsf(p[6], '+ss+', p[1], p[2]) + '+ $
    ;       'p[3] * gausspsf(p[6], '+ss+', p[4], p[5]) + '+ $
    ;       'p[0]*p[7] * gausspsf(p[8], '+ss+', p[1], p[2]) + '+ $
    ;       'p[3]*p[7] * gausspsf(p[8], '+ss+', p[4], p[5])'
endif else if keyword_set(threegauss) then begin
    ; fit using 3-gaussian PSF, using circular gaussians
    if not(keyword_set(silent)) then $
      print, 'using triple-gaussian model for PSF (be extra patient)'
    expr = '     p[0] * gaussian2d('+ss+', p[1], p[2], p[6]/2.35, p[6]/2.35'+ess+', /norm) + '+ $
           'p[3]*p[0] * gaussian2d('+ss+', p[4], p[5], p[6]/2.35, p[6]/2.35'+ess+', /norm) + '+ $
           '     p[0]*p[7] * gaussian2d('+ss+', p[1], p[2], p[8]/2.35, p[8]/2.35'+ess+', /norm) + '+ $
           'p[3]*p[0]*p[7] * gaussian2d('+ss+', p[4], p[5], p[8]/2.35, p[8]/2.35'+ess+', /norm) + '+ $
           '     p[0]*p[9] * gaussian2d('+ss+', p[1], p[2], p[10]/2.35, p[10]/2.35'+ess+', /norm) + '+ $
           'p[3]*p[0]*p[9] * gaussian2d('+ss+', p[4], p[5], p[10]/2.35, p[10]/2.35'+ess+', /norm)'
    ;expr = 'p[0] * gausspsf(p[6], '+ss+', p[1], p[2]) + ' + $
    ;       'p[3] * gausspsf(p[6], '+ss+', p[4], p[5]) + ' + $
    ;       'p[0]*p[7] * gausspsf(p[8], '+ss+', p[1], p[2]) + ' + $
    ;       'p[3]*p[7] * gausspsf(p[8], '+ss+', p[4], p[5]) + ' + $
    ;       'p[0]*p[9] * gausspsf(p[10], '+ss+', p[1], p[2]) + ' + $
    ;       'p[3]*p[9] * gausspsf(p[10], '+ss+', p[4], p[5])'
endif else begin
    ; DEFAULT: fit using circular gaussian PSF
    if not(keyword_set(silent)) then $
      print, 'using single-gaussian model for PSF'
    expr = '     p[0] * gaussian2d('+ss+', p[1], p[2], p[6]/2.35, p[6]/2.35'+ess+', /norm) + '+ $
           'p[3]*p[0] * gaussian2d('+ss+', p[4], p[5], p[6]/2.35, p[6]/2.35'+ess+', /norm)'
    ;expr = 'p[0] * gausspsf(p[6], '+ss+', p[1], p[2]) + '+ $
    ;       'p[3] * gausspsf(p[6], '+ss+', p[4], p[5])'
endelse
if not(keyword_set(silent)) then $
   print, estr


; if desired, add variable for sky level fit
if (keyword_set(skyfit)) then begin
   if not(keyword_set(silent)) then $
      print, 'and fitting for a constant sky level'
   start = [start, 0.0]
   np = n_elements(start)
   expr = expr + ' + p['+strc(n_elements(start)-1)+']'
endif
if (keyword_set(skyslope)) then begin
   if not(keyword_set(silent)) then $
      print, 'and fitting for a sky surface (const + linear gradient)'
   start = [start, 0.0, 0.0, 0.0]
   np = n_elements(start)
   expr = expr + ' + flatsurf2d('+ss+', p['+strc(n_elements(start)-3)+'], '+$
                                       'p['+strc(n_elements(start)-2)+'], '+$
                                       'p['+strc(n_elements(start)-1)+'] )'
endif
if not(keyword_set(skyfit)) and not(keyword_set(skyslope)) and $
   not(keyword_set(silent)) then $
      print, 'not fitting for sky level, assumed to be 0.0'


; fix the (x,y) positions if desired
; --> not quite right, doesn't work for trimmed images!
nparms = n_elements(start)
parinfo=replicate({fixed:0, limited:[0,0], limits:[0.d,0.d]}, nparms)
;parinfo = replicate({fixed:0}, nparms)
if keyword_set(xyfixed) then begin
   start(1:2) = xy1
   start(4:5) = xy2
   parinfo([1, 2, 4, 5]).fixed = 1
endif 

; impose some other common sense limits on the model parameters
; 1. integrated fluxes cannot be negative
parinfo[0].limited(0)=1 & parinfo[0].limits(0)=0.d
;parinfo[3].limited(0)=1 & parinfo[3].limits(0)=0.d
if keyword_set(twogauss) OR keyword_set(threegauss) then begin
   parinfo[7].limited(0)=1 & parinfo[7].limits(0)=0.d
endif
if keyword_set(threegauss) then begin
   parinfo[9].limited(0)=1 & parinfo[9].limits(0)=0.d
endif
; 2. the positions cannot go outside the box
parinfo[1].limited(0)=1 & parinfo[1].limits(0)=0.d
parinfo[2].limited(0)=1 & parinfo[2].limits(0)=0.d
parinfo[4].limited(0)=1 & parinfo[4].limits(0)=0.d
parinfo[5].limited(0)=1 & parinfo[5].limits(0)=0.d
if keyword_set(NOZOOM) then begin
   parinfo[1].limited(1)=1 & parinfo[1].limits(1)=sz[1]-1
   parinfo[2].limited(1)=1 & parinfo[2].limits(1)=sz[2]-1
   parinfo[4].limited(1)=1 & parinfo[4].limits(1)=sz[1]-1
   parinfo[5].limited(1)=1 & parinfo[5].limits(1)=sz[2]-1
endif else begin
   parinfo[1].limited(1)=1 & parinfo[1].limits(1)=zoom-1
   parinfo[2].limited(1)=1 & parinfo[2].limits(1)=zoom-1
   parinfo[4].limited(1)=1 & parinfo[4].limits(1)=zoom-1
   parinfo[5].limited(1)=1 & parinfo[5].limits(1)=zoom-1
endelse
; 3. no FWHM can be less than 1 pixel
parinfo[6].limited(0)=1 & parinfo[6].limits(0)=1.d
if keyword_set(twogauss) OR keyword_set(threegauss) then begin
   parinfo[ 8].limited(0)=1 & parinfo[ 8].limits(0)=1.d
endif
if keyword_set(threegauss) then begin
   parinfo[10].limited(0)=1 & parinfo[10].limits(0)=1.d
endif
; 4. flux ratio of the two components must be 1e-2 < df < 1
parinfo[3].limited(0)=1 & parinfo[3].limits(0)=0.d;1d-4 
parinfo[3].limited(1)=1 & parinfo[3].limits(1)=1d4 

;print,expr
; do the fit & generate resulting model
if not keyword_set(mask) then mask=cut*0+1
if not(keyword_set(nozoom)) then begin
   if keyword_set(xyfixed) then $
      message, '** /xyfixed does not work unless /nozoom is set! **'
   pp = mpfitexpr(expr, indgen(ZOOM, ZOOM), cut, 1.0, start, $
                  weights=mask, parinfo=parinfo, /quiet)
   if total(finite(pp) ne 0) then $
      model = mpevalexpr(expr, indgen(ZOOM, ZOOM), pp) $
   else begin
      message, '**fitting failed**', /info
      pp = start-start
      return
   endelse
endif else begin
   sz = size(cut)
   pp = mpfitexpr(expr, indgen(sz(1), sz(2)), cut, 1.0, start, $
                  weights=mask, parinfo = parinfo, /quiet)
   if total(finite(pp) ne 0) then begin
      model = mpevalexpr(expr, indgen(sz(1), sz(2)), pp) 
   endif else begin 
      message, '**fitting failed**', /info
      pp = start-start
      sep = 0
      pa = 0
      dmag = 0
      return
   endelse
endelse


; define primary as brighter component
; (necessary for BINPARMS.PRO)
if (pp[3] gt 1.0) and not(keyword_set(noswap)) then begin
    message, 'swapping primary and secondary', /info
    tmp = pp[[1, 2]]
    pp[[1, 2]] =  pp[[4, 5]]
    pp[[4, 5]] = tmp
    tmp=1./pp[3]
    pp[0]=pp[3]*pp[0]
    pp[3]=tmp
endif


;------------------------------------------------------------
; show results
;------------------------------------------------------------
;binparms, [pp(1), pp(2)], [pp(4), pp(5)], sep, pa, /silent
sep=      sqrt( (pp[4]-pp[1])^2 + (pp[5]-pp[2])^2 )
pa =posang_pix(  pp[4]-pp[1],      pp[5]-pp[2]    )
dmag = -mag(1./pp(3))
if not(keyword_set(silent)) then begin
    print, '----------------------------------------'
    print, 'RESULTS'
    print, '----------------------------------------'
    print, 'flux ratio (linear units):'
    ;print, '  initial guess = ', strc(start(0)/start(3)), $
    ;       ' ('+strc(mag(start(3)/start(0)))+' mags)'
    print, '  final fit = ', strc(1./pp(3)), $
           ' ('+strc(mag(pp(3)))+' mags)'
    if not(keyword_set(twogauss)) and not(keyword_set(elliptical)) then begin
        print, 'FWHM = ', strc(pp(6)), ' pix'
    endif else begin
        print, 'for double-gaussian PSF:'
        print, '  FWHM(a) = ', strc(pp(6))
        print, '  FWHM(b) = ', strc(pp(8))
        print, '  integrated flux of gauss(a)/gauss(b) = ', 1./pp(7)
        peaka = pp(0) / ((2*!pi) * (pp(6)/2.35)^2.)
        peakb = pp(0) * pp(7) / ((2*!pi) * (pp(7)/2.35)^2.)
        print, '  peak flux of gauss(a)/gauss(b) = ', peaka/peakb
    endelse
    print, 'astrometry:'
    print, '  separation = ', strc(sep), ' pix'
    print, '  position angle = ', strc(pa), ' degs'
    if keyword_set(elliptical) then $
      print, 'elliptical PSF: ratio of FWHMs = ', strc(pp(ii_eratio))
endif

;------------------------------------------------------------
; make rotated versions for plotting crosscut
; rotate about the primary position, for convenience
; (should do this about the photocenter of binary, but too annoying for generating cuts)
if finite(pa) ne 1 then pa=0.
rotpa = pa
xrot = pp(1)
yrot = pp(2)
rotcut = rot(cut, rotpa, 1.0, xrot, yrot, cubic = -0.5, /pivot, miss = 0)
rotmodel = rot(model, rotpa, 1.0, xrot, yrot, cubic = -0.5, /pivot, miss = 0)
;; -- old code for rotating, prior to change in convention for BINPARMS.PRO --
;aaa0 = [pa, 180-abs(pa)]
;aaa = [abs(pa), 180-abs(pa)]
;rotpa = min(aaa, wmin)
;rotpa = sign(aaa0(wmin)) * rotpa
;rotcut = rot(cut, rotpa, 1.0, pp(1), pp(2), cubic = -0.5, /pivot, miss = 0)
;rotmodel = rot(model, rotpa, 1.0, pp(1), pp(2), cubic = -0.5, /pivot, miss = 0)


;------------------------------------------------------------
; display results of image fitting
;   don't open new windows if they are already open
;   (handy when doing multiple calls to BINFIT.PRO)
;------------------------------------------------------------
if not(keyword_set(nodisplay)) then begin

    ;------------------------------------------------------------
    ; display extracted region
    ;------------------------------------------------------------
    device, window_state=winstate
    ; (0) data w/true contours
    if (!d.name eq 'X') then $
      if (winstate(0) eq 0) then win, 0 $
        else wset, 0
    loadct, 15, /silent
    display2, asinh_stretch(cut*mask), tit = 'data w/log contours (asinh stretch)', /silent
    ;display2, cut, tv=tv, tit = 'data w/log contours', /silent
    contour, cut/max(cut), levels=reverse(0.9/2.^findgen(5)), /over, col=0

    ; (1) data w/model contours
    ; overplotting fitted positions on data
    if (!d.name eq 'X') then $
      if (winstate(1) eq 0) then win, 1 $
        else wset, 1
    display2, asinh_stretch(cut*mask), /tv, tit = 'data w/model (log) contours (asinh stretch)', /silent
    ;display2, cut, tv=tv, tit = 'data w/model (log) contours', /silent
    contour, model/max(model), levels=reverse(0.9/2.^findgen(5)), /over, col=0
    oplot, [pp(1), pp(4)], [pp(2), pp(5)], ps = 1, sym = 5, col=0, thick = 3

    ; (2) data - model, with model center positions marked
    if (!d.name eq 'X') then $
      if (winstate(2) eq 0) then win, 2 $
      else wset, 2

    display2, asinh_stretch((cut-model)*mask), /tv, tit = 'data - model (asinh stretch), w/model centers', /silent
    ;display2, cut-model, tv=tv, tit = 'data - model (linear stretch)', /silent
    ;contour, model/max(model), levels=reverse(0.9/2.^findgen(5)), /over, col=0
    oplot, [pp(1), pp(4)], [pp(2), pp(5)], ps = 1, sym = 5, col=0, thick = 3

    ; (8 & 9) debugging check for rotation
    win, 8, /check
    display2, rotcut, /tv, tit = 'rotated data for axis cuts'
    win, 9, /check
    display2, rotmodel, /tv, tit = 'rotated model for axis cuts'
    ;print, rotpa, pp(1), pp(2), pp(4), pp(5)

    ; (3) contours of data vs model
    if (!d.name eq 'X') then $
      if (winstate(3) eq 0) then win, 3 $
        else wset, 3
    display2, cut-cut, /tv, tit = 'data (white) vs model (green)', /silent
    lincolr, /silent
    contour, cut/max(cut), levels=reverse(0.9/2.^findgen(5)), /over, col=3
    contour, model/max(model), levels=reverse(0.9/2.^findgen(5)), /over, col=2

    ;------------------------------------------------------------
    ; plot crosscuts of data vs fitted model
    ;------------------------------------------------------------
    ; (4) major axis cut, linear scale
    if (!d.name eq 'X') then $
      if (winstate(4) eq 0) then win, 4 $
        else wset, 4
    plot, rotcut(xrot, *), ps = 10, tit = 'major axis cut (linear scale)', /xs
    oplot, rotmodel(xrot, *), ps = 10, col = 2
    legend, ['data', 'fit'], col = [!p.color, 2], box = 0, psym = [0, 0]

    ; (5) major axis cut, log scale
    if (!d.name eq 'X') then $
      if (winstate(5) eq 0) then win, 5$
        else wset, 5
    plot, rotcut(xrot, *), ps = 10, tit = 'major axis cut (log scale)', $
          /xs, /ylog, yr=limits(rotcut(xrot, *)) > 1e-2, /ys
    oplot, rotmodel(xrot, *), ps = 10, col = 2
    legend, ['data', 'fit'], col = [!p.color, 2], box = 0, psym = [0, 0]

    ; (6) minor axis cut, linear scale
    if (!d.name eq 'X') then $
      if (winstate(6) eq 0) then win, 6 $
        else wset, 6
    plot, rotcut(*, yrot), ps = 10, tit = 'minor axis cut (linear scale)', $
          yr=limits(rotcut(*, yrot)) > 1e-2, /xs
    oplot, rotmodel(*, yrot), ps = 10, col = 3
    legend, ['data', 'fit'], col = [!p.color, 3], box = 0, psym = [0, 0]

    ; (7) minor axis cut, linear scale
    if (!d.name eq 'X') then $
      if (winstate(7) eq 0) then win, 7 $
        else wset, 7
    plot, rotcut(*, yrot), ps = 10, tit = 'minor axis cut (log scale)', $
          /ylog, yr=limits(rotcut(*, yrot)) > 1e-2, /xs
    oplot, rotmodel(*, yrot), ps = 10, col = 3
    legend, ['data', 'fit'], col = [!p.color, 3], box = 0, psym = [0, 0]

endif


; add offset to return fitted positions to native detector coordinates
pp[[1,4]]=pp[[1,4]]+xoff
pp[[2,5]]=pp[[2,5]]+yoff
; pass fitting results to output variable
fittedparms = pp


end



