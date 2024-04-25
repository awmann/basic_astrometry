
dir = 'data_binfit/'

extr = ''
obj = object+'AB'

; main file containing binfit results
xyfile = file_search( dir+obj+'_????-??-??_*_*'+extr+'.txt', count=nep )
qq = where(strpos(xyfile,'#') eq -1 and strpos(xyfile,'NIRI') eq -1 and strpos(xyfile,'NACO') eq -1)
if qq[0] ne -1 then begin
   xyfile = xyfile[qq]
   nep = n_elements(xyfile)
endif
print,strc(nep)+' binfit epochs found'
forprint,xyfile
print
; (optional) file containing specialized effective wavelengths for the
; two components
dcrfile = file_search( dir+obj+'_effwave.txt', count=ndcr )


mjd_ = -1d
 x0_ = -1.
 y0_ = -1.
 x1_ = -1.
 y1_ = -1.
 df_ = -1.
fil_ = ''
cam_ = ''
ang_ = -1.
wav_ = -1d
pi_ = ''
ut_ = ''
for i=0,nep-1 do begin

   ; read in binfit data from this epoch
   readcol, xyfile[i], f='(d,f,f,f,f,f,a,a,a,a)', $
            mjd, x0, y0, x1, y1, df, filter, camera, imfile,/silent
   ;for k=0,n_elements(imfile)-1 do imfile[k] = strreplace(imfile[k],'gfk','slave1')
   ;for k=0,n_elements(imfile)-1 do imfile[k] = strreplace(imfile[k],'galileo7','slave1')
   ;for k=0,n_elements(imfile)-1 do imfile[k] = strreplace(imfile[k],'slave1','gfk')
   for k=0,n_elements(imfile)-1 do imfile[k] = strreplace(imfile[k],'/Volumes/Vali/','/Volumes/Baldur/')
   for k=0,n_elements(imfile)-1 do imfile[k] = strreplace(imfile[k],'Vali_backup','/Volumes/Baldur/')
   for k=0,n_elements(imfile)-1 do imfile[k] = strreplace(imfile[k],'/Volumes/Loki//Volumes/Loki/','/Volumes/Baldur/')
   ;;for k = 0,n_elements(imfile)-1 do imfile[k] = strreplace(imfile[k],'/Volumes/Loki/','/Volumes/Vali/')
   
   ; toss images flagged as bad
   wbad = where( x0 eq x1 and y0 eq y1, nbad, ncomp=nfile )
   if nbad gt 0 then begin
      remove, wbad, mjd, x0, y0, x1, y1, df
      remove, wbad, filter, camera, imfile
   endif

   ; read in extra header info
   instpa = fltarr(nfile)
   effwav = fltarr(nfile)
   progpi = strarr(nfile)
   utdate = strarr(nfile)
   tstrt = dblarr(nfile)
   tstop = dblarr(nfile)
   thead = dblarr(nfile)
   tmp = fltarr(nfile)
   targra = dblarr(nfile)
   targde = dblarr(nfile)
   for j=0, nfile-1 do begin
      ;;if strpos(imfile[j],'Kraus_nirc2') ne -1 then begin
      ;;   imfile[j] = strreplace(imfile[j],'/Volumes/Loki/Kraus_nirc2/Reduced/
         
      print,imfile[j]
      if file_test(imfile[j]) then begin
         hd = headfits(imfile[j])
         
         if strpos(sxpar(hd,'PARANG'),'#') eq -1 and strpos(sxpar(hd,'ROTPPOSN'),'#') eq -1 and strpos(sxpar(hd,'EL'),'#') eq -1 and strpos(sxpar(hd,'INSTANGL'),'#') eq -1 then $
            instpa[j] = sxpar(hd,'PARANG') + sxpar(hd,'ROTPPOSN') - $
                        sxpar(hd,'EL') - sxpar(hd,'INSTANGL')
                                ;tmp[j] = nirc2derot(imfile[j])
         
         effwav[j] = sxpar(hd,'EFFWAVE')
         
         thead[j] = secstr(sxpar(hd,'UTC'))
         if sxpar(hd,'EXPSTART') eq 0 then begin
            tstrt[j] = thead[j]
            tstop[j] = thead[j]
         endif else begin
            tstrt[j] = secstr(sxpar(hd,'EXPSTART'))
            tstop[j] = secstr(sxpar(hd,'EXPSTOP'))
         endelse
         
                                ;print, mjd[j], '   '+sxpar(hd,'MJD-OBS')
         
         if strc(sxpar(hd,'PROGPI')) ne '0' then $
            progpi[j] = strc(sxpar(hd,'PROGPI')) $
         else $
            progpi[j] = ''
         
         utdate[j] = strc(sxpar(hd,'DATE-OBS'))
         
         targra[j] = sxpar(hd,'RA')
         targde[j] = sxpar(hd,'DEC')
         
                                ; clean up filter names
         filter[j] = strreplace(strreplace(strreplace(strreplace(filter[j],'PK50_1.5',''),'clear',''),'_',''),'+','')
      endif
   endfor
      
   if nfile gt 1 then $
      instpa2 = interpol( instpa, thead, (tstop+tstrt)/2d ) $
   else $
      instpa2 = instpa
   ;print,xyfile[i]
   ;print,instpa2-instpa
   ;print,tmp
   ;stat,instpa2-instpa-tmp
   ;print,separation(x1-x0,y1-y0)
   ;print,posang_pix(x1-x0,y1-y0)+instpa2+360
   ;stop

   mjd_ = [ mjd_, mjd ]
    x0_ = [  x0_,  x0 ]
    y0_ = [  y0_,  y0 ]
    x1_ = [  x1_,  x1 ]
    y1_ = [  y1_,  y1 ]
    df_ = [  df_,  df ]
   fil_ = [ fil_, filter ]
   cam_ = [ cam_, camera ]
   ang_ = [ ang_, instpa2 ]
   wav_ = [ wav_, effwav ]
    pi_ = [  pi_, progpi ]
    ut_ = [  ut_, utdate ]

endfor
mjd_ = mjd_[1:*]
 x0_ =  x0_[1:*]
 y0_ =  y0_[1:*]
 x1_ =  x1_[1:*]
 y1_ =  y1_[1:*]
 df_ =  df_[1:*]
fil_ = fil_[1:*]
cam_ = cam_[1:*]
ang_ = ang_[1:*]
wav_ = wav_[1:*]
 pi_ =  pi_[1:*]
 ut_ =  ut_[1:*]
nmeas = n_elements(mjd_)

se_ins = fltarr(nmeas)
pa_ins = fltarr(nmeas)
se_sky = fltarr(nmeas)
pa_sky = fltarr(nmeas)
w = where( cam_ eq 'narrow', n )
if n gt 0 then begin

   print, 'applying narrow camera distortion solution...'

   ; apply distortion correction
   nirc2_dedistort_xy, x0_[w], y0_[w], x0c, y0c, camera='narrow'
   nirc2_dedistort_xy, x1_[w], y1_[w], x1c, y1c, camera='narrow', $
                       scale=scl, orient=rot, escale=escl, eorient=erot

   ; check if there is data from after the AO realignment
   postrealign = 0
   w2 = where( mjd_[w] gt 57125d, n2 )
   if n2 gt 0 then begin
      nirc2_dedistort_xy, x0_[w[w2]], y0_[w[w2]], x0c2, y0c2, camera='narrow', /postrealign
      nirc2_dedistort_xy, x1_[w[w2]], y1_[w[w2]], x1c2, y1c2, camera='narrow', /postrealign, scale=scl2, orient=rot2, escale=escl2, eorient=erot2
      x0c[w[w2]] = x0c2
      y0c[w[w2]] = y0c2
      x1c[w[w2]] = x1c2
      y1c[w[w2]] = y1c2
      postrealign = 1      
   endif

   ; compute instrument sep and PA
   se_ins[w] = separation( x1c-x0c, y1c-y0c )
   pa_ins[w] = posang_pix( x1c-x0c, y1c-y0c )

   ; compute sky sep and PA
   se_sky[w] = se_ins[w] * scl
   pa_sky[w] = pa_ins[w] + rot + ang_[w]
   if postrealign eq 1 then begin
      se_sky[w[w2]] = se_ins[w[w2]] * scl2
      pa_sky[w[w2]] = pa_ins[w[w2]] + rot2 + ang_[w[w2]]
   endif

endif
;forprint,mjd_[w],x0c,y0c,x1c,y1c,ang_[w] & stop
;
w = where( cam_ eq 'wide', n )
if n gt 0 then begin

   print, 'applying wide camera distortion solution...'

   ; apply distortion correction
   nirc2_dedistort_xy, x0_[w], y0_[w], x0c, y0c, camera='wide'
   nirc2_dedistort_xy, x1_[w], y1_[w], x1c, y1c, camera='wide', $
                       scale=scl, orient=rot

   ; compute instrument sep and PA
   se_ins[w] = separation( x1c-x0c, y1c-y0c )
   pa_ins[w] = posang_pix( x1c-x0c, y1c-y0c )

   ; compute sky sep and PA
   se_sky[w] = se_ins[w] * scl
   pa_sky[w] = pa_ins[w] + rot + ang_[w]

endif

; compute Delta RA and Dec between components
dr_sky = se_sky * sin(pa_sky/!radeg)
dd_sky = se_sky * cos(pa_sky/!radeg)

; determine atmospheric conditions
print, 'finding weather data...'
if n_elements(temp_c) ne 3381653l then $
   restore, '~/Dropbox/BINFIT/mko-weather.sav'
;w = where( mjd_mko ge min(mjd_) and mjd_mko le max(mjd_) )
;mjd_mko  =  mjd_mko[w]  
;wspd_kts = wspd_kts[w] 
;wdir_deg = wdir_deg[w] 
;temp_c   =   temp_c[w]   
;relh_pct = relh_pct[w] 
;pres_mb  =  pres_mb[w]  
temp = fltarr(nmeas)
pres = fltarr(nmeas)
relh = fltarr(nmeas)
for i=0,nmeas-1 do begin
   w1 = where(finite(temp_c))    &  tmp1 = min( abs(mjd_mko[w1]-mjd_[i]), w1b )
   w2 = where(finite(pres_mb))   &  tmp2 = min( abs(mjd_mko[w2]-mjd_[i]), w2b )
   w3 = where(finite(relh_pct))  &  tmp3 = min( abs(mjd_mko[w3]-mjd_[i]), w3b )
   temp[i] = double(temp_c[w1[w1b]])
   pres[i] = double(pres_mb[w2[w2b]]*100d)
   relh[i] = double(relh_pct[w3[w3b]])
   if tmp1 gt 0.5/24. or tmp2 gt 0.5/24. or tmp3 gt 0.5/24. then $
      message,string([tmp1,tmp2,tmp3]*24.,f='(3(f5.1,x))')+string(mjd_[i],format="(D20.5)")+string(9b)+date[i],/cont
endfor
w = where( finite(pres) ne 1, comp=w2, n )
if n gt 0 then pres[w] = median(pres[w2])
w = where( finite(temp) ne 1, comp=w2, n )
if n gt 0 then temp[w] = median(relh[w2])
w = where( finite(relh) ne 1, comp=w2, n )
if n gt 0 then relh[w] = median(relh[w2])


; compute DAR offsets
print, 'computing DAR...'
wav1_ = wav_
wav2_ = wav_
;
if ndcr eq 1 then begin
   readcol, dcrfile[0], tmpfilt, tmpl1, tmpl2, f='(a,f,f)',/silent
   for i=0,n_elements(tmpfilt)-1 do begin
      wfil = where( fil_ eq tmpfilt[i], ntmp )
      if ntmp gt 0 then begin
         print,'using different effective wavelengths for primary and secondary in '+tmpfilt[i]
         wav1_[wfil] = tmpl1[i]
         wav2_[wfil] = tmpl2[i]
      endif
   endfor
endif
;
refraction, 2400000.5d + mjd_, $
            median(targra)+dblarr(nmeas), $
            median(targde)+dblarr(nmeas), $
            temp, pres, relh, wav1_, 'keck', ra_refr_off0, de_refr_off0, zen=zen_
refraction, 2400000.5d + mjd_, $
            median(targra)+dr_sky/3600e3/cos(median(targde)/!radeg), $
            median(targde)+dd_sky/3600e3, $
            temp, pres, relh, wav2_, 'keck', ra_refr_off1, de_refr_off1


; compute aberration offsets
print, 'computing aberration...'
aberration, 2400000.5d + mjd_, $
            median(targra)+dblarr(nmeas), $
            median(targde)+dblarr(nmeas), $
            ra_aber_off0, de_aber_off0
aberration, 2400000.5d + mjd_, $
            median(targra)+dr_sky/3600e3/cos(median(targde)/!radeg), $
            median(targde)+dd_sky/3600e3, $
            ra_aber_off1, de_aber_off1


; compute final sep and PA
tmpdr = dr_sky-(((ra_refr_off1-ra_refr_off0)+(ra_aber_off1-ra_aber_off0))*3600e3/cos(median(targde)/!radeg))
tmpdd = dd_sky-((de_refr_off1-de_refr_off0)+(de_aber_off1-de_aber_off0))*3600e3
se_fin = separation( -tmpdr, tmpdd )
pa_fin = posang_pix( -tmpdr, tmpdd )


; deal with 0/360 boundary
; pa_fin = dphase(pa_fin,max(pa_fin,/nan))+max(pa_fin,/nan)
;forprint,pa_ins,pa_sky+360,pa_fin

; find unique epochs and save results for each
epmjd = long(mjd_)
uepmjd = epmjd[uniq(epmjd,sort(epmjd))]
nuep = n_elements(uepmjd)
dravg = -1.
ddavg = -1.
seavg = -1.
paavg = -1.
dmavg = -1.
drrms = -1.
ddrms = -1.
serms = -1.
parms = -1.
dmrms = -1.
spcov = -1.
smcov = -1.
pmcov = -1.
rdcov = -1.
rmcov = -1.
dmcov = -1.
mjd_u = -1d
epstr = '' 
pistr = '' 
fistr = '' 
ndith = -1 
for i=0,nuep-1 do begin

   wep = where( epmjd eq uepmjd[i] )
   tmp = fil_[wep] 
   ufilt = tmp[uniq(tmp,sort(tmp))]
   nfilt = n_elements(ufilt)

   for j=0,nfilt-1 do begin

      wf = where( fil_[wep] eq ufilt[j], nstart )

      ; kludge to handle single measurements
      if nstart eq 1 then begin
         wf = [ wf[0], wf[0] ]
         nstart = 2
      endif

      if nstart gt 1 then begin

      tmps = se_fin[wep[wf]]
      tmpp = pa_fin[wep[wf]]
      ;tmpf = df_[wep[wf]]
      tmpf = 2.5*alog10(df_[wep[wf]])
      tmpi = indgen(nstart)
      tmpr = tmps*sin(tmpp/!radeg)
      tmpd = tmps*cos(tmpp/!radeg)

      nclip = 999
      ngood = 999
      if nstart lt 20 then $
         sig = 3.0
      if nstart ge 20 then $
         sig = 3.5
      while nclip gt 0 and ngood gt 2 do begin
         w = where( abs(tmpr-median(tmpr,/even)) le sig*stddev(tmpr,/nan) and $
                    abs(tmpd-median(tmpd,/even)) le sig*stddev(tmpd,/nan) and $
                    abs(tmps-median(tmps,/even)) le sig*stddev(tmps,/nan) and $
                    abs(tmpp-median(tmpp,/even)) le sig*stddev(tmpp,/nan) and $
                    abs(tmpf-median(tmpf,/even)) le sig*stddev(tmpf,/nan), $
                    ngood, comp=tmpw, ncomp=nclip )
         ;stop
         tmpr = tmpr[w]
         tmpd = tmpd[w]
         tmps = tmps[w]
         tmpp = tmpp[w]
         tmpf = tmpf[w]
         tmpi = tmpi[w]
      endwhile
      print, 'ep'+strc(i)+' '+string(ufilt[j],f='(a4)')+': clipped '+strc(nstart-ngood)+' of '+strc(nstart)+' measurements'

                                ; flip some systems
      if n_elements(tmpp) gt 1 then begin
         w = where(abs(dphase(tmpp,median(tmpp))) gt 90, nflip)
         if nflip gt 0 then begin
            tmpr[w] = -tmpr[w]
            tmpd[w] = -tmpd[w]
            tmpp[w] = (tmpp[w]+180.) mod 360
                                ;tmpf[w] = 1./tmpf[w]
            tmpf[w] = -tmpf[w]
            print, 'ep'+strc(i)+' '+string(ufilt[j],f='(a4)')+': flipped '+strc(nflip)+' of '+strc(ngood)+' measurements'
         endif
      endif

      ; deal with 0/360 boundary
      tmpp = dphase(tmpp,min(tmpp,/nan))+((min(tmpp,/nan)+360.) mod 360.)

      ; flip special cases
      if (obj eq 'SD0926+58' and fil_[wep[wf[tmpi[0]]]] ne 'J') or $
         (obj eq '2M1404-31' and fil_[wep[wf[tmpi[0]]]] eq 'J') or $
         (obj eq 'SD1021-03' and fil_[wep[wf[tmpi[0]]]] eq 'J') or $
         (obj eq 'SD1534+16' and fil_[wep[wf[tmpi[0]]]] eq 'J') or $
         (obj eq '2M1404-31' and fil_[wep[wf[tmpi[0]]]] eq 'K' and ut_[wep[wf[tmpi[0]]]] eq '2013-04-29') or $
         (obj eq 'DE1228-15' and fil_[wep[wf[tmpi[0]]]] eq 'Y' and ut_[wep[wf[tmpi[0]]]] eq '2013-01-17') or $
         (obj eq 'SD2052-16' and fil_[wep[wf[tmpi[0]]]] eq 'J' and ut_[wep[wf[tmpi[0]]]] eq '2005-10-11') or $
         (obj eq 'SD1052+44' and fil_[wep[wf[tmpi[0]]]] eq 'H' and ut_[wep[wf[tmpi[0]]]] eq '2008-11-03') or $
         (obj eq '2M1743+58' and fil_[wep[wf[tmpi[0]]]] eq 'Y' and ut_[wep[wf[tmpi[0]]]] eq '2013-07-01') or $
         (obj eq 'SD1052+44' and fil_[wep[wf[tmpi[0]]]] eq 'J') then begin
         tmpr = -tmpr
         tmpd = -tmpd
         tmpp = (tmpp+180.) mod 360
         tmpf = -tmpf
         print,'flipping all measurements in this filter/epoch!'
      endif

      ; calculate mean, rms, covariance, etc.
      dravg = [ dravg, mean(tmpr,/nan) ]
      ddavg = [ ddavg, mean(tmpd,/nan) ]
      seavg = [ seavg, mean(sqrt(tmpr^2+tmpd^2),/nan) ]
      paavg = [ paavg, mean(tmpp,/nan) ]
      ;paavg = [ paavg, mean((atan(tmpr,tmpd)*!radeg+360.) mod 360.,/nan) ]
      ;seavg = [ seavg, mean(tmps,/nan) ]
      ;paavg = [ paavg, mean(tmpp,/nan) ]
      dmavg = [ dmavg, mean(tmpf,/nan) ]
      tmpxy = stddev([tmpr-mean(tmpr),tmpd-mean(tmpd)],/nan)
      drrms = [ drrms, tmpxy ]
      ddrms = [ ddrms, tmpxy ]
      ;serms = [ serms, stddev(sqrt(tmpr^2+tmpd^2),/nan) ]
      ;parms = [ parms, stddev((atan(tmpr,tmpd)*!radeg+360.) mod 360.,/nan) ]
      ;drrms = [ drrms, stddev(tmpr,/nan) ]
      ;ddrms = [ ddrms, stddev(tmpd,/nan) ]
      serms = [ serms, stddev(tmps,/nan) ]
      parms = [ parms, stddev(tmpp,/nan) ]
      dmrms = [ dmrms, stddev(tmpf,/nan) ]
      spcov = [ spcov, mean((tmps-mean(tmps,/nan))*(tmpp-mean(tmpp,/nan)),/nan) ]
      smcov = [ smcov, mean((tmps-mean(tmps,/nan))*(tmpf-mean(tmpf,/nan)),/nan) ]
      pmcov = [ pmcov, mean((tmpp-mean(tmpp,/nan))*(tmpf-mean(tmpf,/nan)),/nan) ]
      rdcov = [ rdcov, mean((tmpr-mean(tmpr,/nan))*(tmpd-mean(tmpd,/nan)),/nan) ]
      rmcov = [ rmcov, mean((tmpr-mean(tmpr,/nan))*(tmpf-mean(tmpf,/nan)),/nan) ]
      dmcov = [ dmcov, mean((tmpd-mean(tmpd,/nan))*(tmpf-mean(tmpf,/nan)),/nan) ]
      ;
      mjd_u = [ mjd_u, mean(mjd_[wep[wf[tmpi]]],/nan) ]
      epstr = [ epstr,  ut_[wep[wf[tmpi[0]]]] ]
      pistr = [ pistr,  pi_[wep[wf[tmpi[0]]]] ]
      fistr = [ fistr, fil_[wep[wf[tmpi[0]]]] ]
      ndith = [ ndith, total(finite(tmps)) ]
      ;stop

      endif


   endfor

   if nfilt gt 1 then begin
      print
      print, wgtavg(dravg[n_elements(serms)-nfilt:*],drrms[n_elements(serms)-nfilt:*],chi2=chi2,pchi2=pchi2,adderr=adderr,wavgerr=wavgerr), $
             wavgerr, adderr, pchi2, f='("delr: ",f8.3,x,f7.3,"   adderr=",f7.3,"    pchi2=",f6.4)'
      print, wgtavg(ddavg[n_elements(serms)-nfilt:*],ddrms[n_elements(serms)-nfilt:*],chi2=chi2,pchi2=pchi2,adderr=adderr,wavgerr=wavgerr), $
             wavgerr, adderr, pchi2, f='("deld: ",f8.3,x,f7.3,"   adderr=",f7.3,"    pchi2=",f6.4)'
      print, wgtavg(seavg[n_elements(serms)-nfilt:*],serms[n_elements(serms)-nfilt:*],chi2=chi2,pchi2=pchi2,adderr=adderr,wavgerr=wavgerr), $
             wavgerr, adderr, pchi2, f='("sep.: ",f8.3,x,f7.3,"   adderr=",f7.3,"    pchi2=",f6.4)'
      print, wgtavg(paavg[n_elements(serms)-nfilt:*],parms[n_elements(serms)-nfilt:*],chi2=chi2,pchi2=pchi2,adderr=adderr,wavgerr=wavgerr), $
             wavgerr, adderr, pchi2, f='("p.a.: ",f8.3,x,f7.3,"   adderr=",f7.3,"    pchi2=",f6.4)'
      print, wgtavg(dmavg[n_elements(serms)-nfilt:*],dmrms[n_elements(serms)-nfilt:*],chi2=chi2,pchi2=pchi2,adderr=adderr,wavgerr=wavgerr), $
             wavgerr, adderr, pchi2, f='("dmag: ",f8.3,x,f7.3,"   adderr=",f7.3,"    pchi2=",f6.4)'
      print
   endif

   ;stop

endfor
;stop
dravg = dravg[1:*]
ddavg = ddavg[1:*]
seavg = seavg[1:*]
paavg = paavg[1:*]
dmavg = dmavg[1:*]
drrms = drrms[1:*]
ddrms = ddrms[1:*]
serms = serms[1:*]
parms = parms[1:*]
dmrms = dmrms[1:*]
spcov = spcov[1:*]
smcov = smcov[1:*]
pmcov = pmcov[1:*]
rdcov = rdcov[1:*]
rmcov = rmcov[1:*]
dmcov = dmcov[1:*]
mjd_u = mjd_u[1:*]
epstr = epstr[1:*]
pistr = pistr[1:*]
fistr = fistr[1:*]
ndith = ndith[1:*]

forprint, text='astrom-'+obj+extr+'-sp.txt', f='(a-10,2x,f10.2,4x,f8.3,x,f7.3,4x,f8.3,x,f7.3,3x,f8.4,x,f7.4,3x,i3,3x,a-4,3x,3f8.3,3x,a)', $
          epstr, mjd_u, seavg, serms, paavg, parms, dmavg, dmrms, $
          ndith, fistr, spcov, smcov, pmcov, pistr, subset=sort(mjd_u)

forprint, text='astrom-'+obj+extr+'-xy.txt', f='(a-10,2x,f10.2,4x,f9.3,x,f7.3,4x,f9.3,x,f7.3,3x,f8.4,x,f7.4,3x,i3,3x,a-4,3x,3f8.3,3x,a)', $
          epstr, mjd_u, dravg, drrms, ddavg, ddrms, dmavg, dmrms, $
          ndith, fistr, rdcov, rmcov, dmcov, pistr, subset=sort(mjd_u)

;;if not exist('extra_info/'+obj+extr+'.exclude') then $
forprint, text='extra_info/'+obj+extr+'.exclude', f='(a-12,3x,a-4,3x,a)', $
          '# '+epstr, fistr, pistr, subset=sort(mjd_u)

;;if not exist('extra_info/offsets-'+obj+extr+'-sp.txt') then $
forprint, text='extra_info/offsets-'+obj+extr+'-sp.txt', f='(a-10,2x,f10.2,4x,f8.3,x,f7.3,4x,f8.3,x,f7.3,3x,f8.4,x,f7.4,3x,i3,3x,a-4,3x,a)', $
          epstr, mjd_u, seavg*0., serms*0., paavg*0., parms*0., dmavg*0., dmrms*0., $
          ndith, fistr, pistr, subset=sort(mjd_u)

;;if not exist('extra_info/offsets-'+obj+extr+'-xy.txt') then $
forprint, text='extra_info/offsets-'+obj+extr+'-xy.txt', f='(a-10,2x,f10.2,4x,f8.3,x,f7.3,4x,f8.3,x,f7.3,3x,f8.4,x,f7.4,3x,i3,3x,a-4,3x,a)', $
          epstr, mjd_u, dravg*0., drrms*0., ddavg*0., ddrms*0., dmavg*0., dmrms*0., $
          ndith, fistr, pistr, subset=sort(mjd_u)

end
