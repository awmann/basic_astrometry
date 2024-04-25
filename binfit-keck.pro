

drive = 'Vali'

comp0 = 0
comp1 = 1

do_fit = 1
sf_fit = 0;; starfinder (=1) versus regular (=0), starfinder usually better but higher fail rate?
zoommult = 2.5;; normally ~2.5. This decides how 'big' a window you are using. If you find the PSFs are not 'inside' the X11 windows, increase the size a bit.
;;wide_zoom = 56
wide_zoom = -1 ;; if this is >0 then this is the number of pixels in the window. Use 50-100 or so, depends on target. Only >0 if the targets are widely separated (like >>1" or thereabouts). This is going to be massive trial and error, 50 might be perfect or might crash.
qzap = 0
slope= 0
; for SFBINFIT
postol = 0.010 ;; tolerance for convergence in pixels
flxtol = 0.005 ;; tolerance for convergence in flux or maybe mags, Trent isn't sure, which is why he should comment his code'


if comp0 eq 0 then comp0str='A'
if comp1 eq 1 then comp1str='B'


; find all data files
find_keck_data, object, filename=files, /mask, drive=drive

nfiles = n_elements(files)
drive = 'Vali'
for i=0,nfiles-1 do begin
   files[i] = strreplace(files[i],'/Volumes/Vali/','/Volumes/Loki/')
endfor
;;for i = 0,n_elements(files)-1 do files[i] = strreplace(files[i],'/Volumes/Loki/','/Volumes/Vali/')


; read in blacklist
readcol, 'blacklist.txt', blackfiles, f='(a)', /silent, comment='#'
excl = bytarr(nfiles)
for i=0,nfiles-1 do begin
   fdecomp, files[i], disk__, dir__, name__, qual__
   w = where( strpos(blackfiles,name__) ne -1, n )
   if n gt 0 then excl[i] = 1b
endfor
print, strc(fix(total(excl)))+' of '+strc(nfiles)+' files on blacklist'
print
; remove any files from consideration as needed
w = where( excl, n )
if n gt 0 then $
   remove, w, files
nfiles = n_elements(files)


; identify groupings from header info
cam    = strarr(nfiles)
filt   = strarr(nfiles)
date   = strarr(nfiles)
instpa = fltarr(nfiles)
t0     = dblarr(nfiles)
;
for i=0,nfiles-1 do begin

   hd = headfits( files[i], /silent )

   cam[i]  = strc(sxpar(hd,'CAMNAME'))
   filt[i] = strc(sxpar(hd,'FILTER'))
   date[i] = strc(sxpar(hd,'DATE-OBS'))

   instpa[i] = sxpar(hd,'PARANG') + sxpar(hd,'ROTPPOSN') - $
               sxpar(hd,'EL')     - sxpar(hd,'INSTANGL')

   t0[i] = sxpar(hd,'MJD-OBS')

                                ; do some cosmetic stuff to filter name
   filt[i] = strreplace(filt[i],'PK50_1.5','')
   filt[i] = strreplace(filt[i],'clear','')
   filt[i] = strreplace(filt[i],'_','')
   filt[i] = strreplace(filt[i],'+','')

endfor
;
groupstr = strc(date)+'_'+strc(filt)+'_'+strc(cam)
ugroup = groupstr[uniq(groupstr,sort(groupstr))]
ngroup = n_elements(ugroup)


; figure out which data sets to analyze
for i=0,ngroup-1 do begin
   w = where( groupstr eq ugroup[i], nimg )
   print, i, ugroup[i], nimg, f='(i3,2x,a-30,i3," files")'
endfor
ii = lindgen(ngroup)
if getyn('analyze all non-masking data?') ne 1 then begin
   tmp = 0
   read, prompt='which one to analyze? ', tmp
   ii = tmp
endif


;hacks go here
;ii = [ 3, 4, 5, 6 ] & print,'HACK!!!!!!!!!'


; perform binfit analysis for each group
for k=0,n_elements(ii)-1 do begin

   i = ii[k]

   if strpos(ugroup[i],'hole') eq -1 then begin

                                ; some strings we'll be using'
      base = object+comp0str+comp1str+'_'+ugroup[i]
      qzs = '.qzap' & if qzap ne 1 then qzs=''
      sks = '.sky' & if slope eq 1 then sks='.skyslope.sky'

                                ; find files in this group
      w = where( groupstr eq ugroup[i], nimg )
      print, i, ugroup[i], nimg, f='(i3,2x,a-30,i3," files")'

                                ; load images into data cube
      xadd = dblarr(nimg)
      yadd = xadd
      for j=0,nimg-1 do begin
         tmp = readfits(files[w[j]],/silent)
         if n_elements(tmp[*,0]) lt 1024 or n_elements(tmp[0,*]) lt 1024 then begin
            xadd[j] = (1024- n_elements(tmp[*,0]))/2.
            yadd[j] = (1024- n_elements(tmp[0,*]))/2.
            old = tmp
            tmp = reform_nirc2(tmp)
            tmp[where(tmp le 0)] = median(old)
         endif
         if j eq 0 then imcube = tmp else imcube = [ [[imcube]], [[tmp]] ]
      endfor
      if nimg eq 1 then imcube = reform( imcube, [ size(imcube,/dim), 1 ] )

      ;; find info on xy locations and save as needed
      tmp = strarr(nimg)
      x0 = fltarr(nimg)
      y0 = fltarr(nimg)
      dx = fltarr(nimg)
      dy = fltarr(nimg)
      for j=0,nimg-1 do begin
         fdecomp, files[w[j]], disk__, dir__, name__, qual__
         readcol, 'centers/'+name__+'.xy.txt', tmpx, tmpy, /silent
         x0[j] = tmpx[comp0]
         y0[j] = tmpy[comp0]
         dx[j] = tmpx[comp1] - tmpx[comp0]
         dy[j] = tmpy[comp1] - tmpy[comp0]
      endfor
      wflip = where( sign(dx) ne sign(dx[0]) and sign(dy) ne sign(dy[0]), nflip )
      if nflip gt 0 then begin
         x0[wflip] = x0[wflip] + dx[wflip]
         y0[wflip] = y0[wflip] + dy[wflip]
         dx[wflip] = -dx[wflip]
         dy[wflip] = -dy[wflip]
      endif
      x0+=xadd
      y0+=yadd

      xyfile = 'binfit_data/'+base+'_shifts.good.dat'
      forprint, text=xyfile, x0, y0, f='(2f10.3,x,"   NaN   1")', /silent

      if wide_zoom lt 0 then begin
         ;; NORMAL CASE

         zoom = round(mean(separation(dx,dy)+2.)*zoommult)
         psfsz = round(zoom*1.0) ;0.7
         if sf_fit eq 1 then $
            results_file = 'binfit_data/sfbinfit_'+base+'.RESULTS' $
         else $
            results_file = 'binfit_data/binfit_'+base+'.3gauss'+sks+qzs+'.'+strc(zoom)+'pix'

                                ; run binfit
         wclear
         if do_fit eq 1 and sf_fit eq 1 then $
            sfbinfit_go, imcube, outstr='binfit_data/sfbinfit_'+base, $
                         postol=postol, flxtol=flxtol, $
                         xyfile=xyfile, dxy_sec=[[dx],[dy]], $
                         zoom=zoom,psfsize=psfsz,/nocheck,qzap=qzap,/good
         if do_fit eq 1 and sf_fit ne 1 then $
            binfit_go, imcube, 'binfit_data/binfit_'+base, $
                       xyfile=xyfile, dxy_sec=[[dx],[dy]], $
                                ;/two,/ellip,/skyfit,skyslope=slope,zoom=zoom,qzap=qzap,/good
                       /three,/ellip,/skyfit,skyslope=slope,zoom=zoom,qzap=qzap,/good


      endif else begin
         ;; WIDE CASE

         if sf_fit eq 1 then $
            results_file = 'binfit_data/sfbinfit_'+base+'.RESULTS' $
         else $
            results_file = 'binfit_data/binfit_'+base+'.3gauss'+sks+qzs+'.wide_zoom'+strc(wide_zoom)+'pix'

         psfsz = wide_zoom
         fitimg = fltarr(psfsz*2,psfsz)
         xorig0 = round(x0)-psfsz/2
         yorig0 = round(y0)-psfsz/2
         xorig1 = round(x0+dx)-psfsz/2
         yorig1 = round(y0+dy)-psfsz/2
         im = lindgen(nimg)
         sep = fltarr(nimg)
         pa = fltarr(nimg)
         df = fltarr(nimg)
         x0 = fltarr(nimg)
         y0 = fltarr(nimg)
         x1 = fltarr(nimg)
         y1 = fltarr(nimg)
         f0 = fltarr(nimg)
         f1 = fltarr(nimg)
         niter = intarr(nimg)
         good = intarr(nimg)+1
         parms = fltarr(nimg,14)
         for j=0,nimg-1 do begin
            if xorig0[j] lt 0 then xorig0[j] = 0

            ;; stitch together an image
            y2 = yorig1[j]+psfsz-1
            ;;if yorig1[j]+psfsz-1 gt n_elements(imcube[0,*,0])-1 then y2 = n_elements(imcube[0,*,0])-1
            fitimg[0    :psfsz-1,*] = imcube[xorig0[j]:xorig0[j]+psfsz-1,yorig0[j]:yorig0[j]+psfsz-1,j]
            fitimg[psfsz:*      ,*] = imcube[xorig1[j]:xorig1[j]+psfsz-1,yorig1[j]:y2,j]
            noiseimg = fitimg*0.+1.

                                ; run binfit
            wclear
            if do_fit eq 1 and sf_fit eq 1 then begin
               print,'----------------------------'
               print,'image '+strc(j+1)+' of '+strc(nimg)
               print,'----------------------------'
               sfbinfit, fitimg, noiseimg, /nozoom, $
                         sep_, pa_, fratio_, niter_, $
                         im_x, im_y, im_flux, $
                         psf, synfield, $
                         postol=postol, flxtol=flxtol, $
                         xy1=[psfsz*0.5,psfsz*0.5], $
                         xy2=[psfsz*1.5,psfsz*0.5], $
                         psfsize=psfsz, outstr='tmp-sfbinfit.dat'
               df[j] = fratio_
               x0[j] = xorig0[j] + im_x[0]
               y0[j] = yorig0[j] + im_y[0]
               x1[j] = xorig1[j] + im_x[1] - psfsz
               y1[j] = yorig1[j] + im_y[1]
               niter[j] = niter_
               f0[j] = im_flux[0]
               f1[j] = im_flux[1]
               sep[j] = separation(x1[j]-x0[j],y1[j]-y0[j])
               pa[j] = posang_pix(x1[j]-x0[j],y1[j]-y0[j])
            endif
            if do_fit eq 1 and sf_fit ne 1 then begin
               print,'----------------------------'
               print,'image '+strc(j+1)+' of '+strc(nimg)
               print,'----------------------------'
               binfit, fitimg, sep_, pa_, parms_, model_, dmag_, $
                       /three, /ellip, /skyfit, skyslope=slope, qzap=qzap, $
                       xy1=[psfsz*0.5,psfsz*0.5], xy2=[psfsz*1.5,psfsz*0.5]
               df[j] = 10.^(0.4*dmag_)
               x0[j] = xorig0[j] + parms_[1]
               y0[j] = yorig0[j] + parms_[2]
               x1[j] = xorig1[j] + parms_[4] - psfsz
               y1[j] = yorig1[j] + parms_[5]
               parms[j,*] = parms_
               sep[j] = separation(x1[j]-x0[j],y1[j]-y0[j])
               pa[j] = posang_pix(x1[j]-x0[j],y1[j]-y0[j])
            endif

         endfor

                                ; make traditional output files
         if do_fit eq 1 and sf_fit eq 1 then $
            forprint, text=results_file, comment='# img   sep(pix)   instrPA    f_0/f_1  niter  good       x_0           y_0           x_1           y_1          f_0          f_1', $
                      f='(a,f10.5,f10.3,f11.5,i6,i6,4f14.6,2e13.3)', $
                      strc(im), sep, pa, df, niter, good, x0, y0, x1, y1, f0, f1
         if do_fit eq 1 and sf_fit ne 1 then $
            forprint, text=results_file, comment='#im  sep(pix)   instrPA    f_0/f_1  good    f_0(A)       f_1(A)     f(B)/f(A)    f(C)/f(A)   FWHMx(A)  FWHMx(B)  FWHMx(C)  FWHMy/FWHMx      x_0           y_0           x_1           y_1   ', $
                      f='(a3,f10.5,f10.3,f11.5,i4,4e13.3,3f10.4,f12.6,4f14.6)', $
                      strc(im), sep, pa, df, good, $
                      parms[*,0], parms[*,3]*parms[*,0], $
                      parms[*,7], parms[*,9], parms[*,6], parms[*,8], parms[*,10], parms[*,11], $
                      x0, y0, x1, y1

      endelse

                                ; determine if this group is in a subarray
      xoff = (1024-n_elements(imcube[*,0,0]))/2
      yoff = (1024-n_elements(imcube[0,*,0]))/2


                                ; read in results
      if sf_fit eq 1 then $
         readcol, results_file, /silent, $
                  im, sep, pa, df, niter, flag, x0, y0, x1, y1, f0, f1 $
      else $
         readcol, results_file, /silent, $
                  im, sep, pa, df, flag, f0a, f1a, fba, fca, $
                  fwhmxa, fwhmxb, fwhmxc, fwhmyx, x0, y0, x1, y1


                                ; output final results in common format
      yn = 1
      if file_test('data_binfit/'+base+'.txt') then $
         yn=getyn('overwrite existing binfit results?')
      if yn then $
         forprint, text='data_binfit/'+base+'.txt', comment=tmp, $
                   f='(f11.5,3x,f9.4,x,f9.4,2x,f9.4,x,f9.4,2x,f9.5,2x,a-'+strc(max(strlen(filt[w])))+',2x,a-'+strc(max(strlen(cam[w])))+',2x,a-'+strc(max(strlen(files[w])))+')', $
                   t0[w], x0+xoff, y0+yoff, x1+xoff, y1+yoff, df, filt[w], cam[w], files[w]
      spawn,'cat data_binfit/'+base+'.txt'

      forprint, sep, pa, 2.5*alog10(df), files[w], f='(3f10.3,3x,a)',comment='Sep(") PA(deg) Fa/Fb path'
      print,'Errors:'
      print, stddev(sep), stddev(pa), stddev(2.5*alog10(df)), f='(3f10.3)'
      stat, imcube
      if getyn('continue?') eq 0 then stop

   endif

endfor

close, /all

end
