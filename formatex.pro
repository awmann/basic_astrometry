
PRO comment,arr,date
  l = where(strpos(arr,date) ne -1)
  if l[0] eq -1 then begin
     print,'Not found!'
     stop
  endif else begin
     arr[l] = '#'+arr[l]
  endelse
END



PRO getmnstr,mn,mnstr
  
  mnstr = strarr(n_elements(mn))
  for i = 0,n_elements(mn)-1 do begin
     case mn[i] of
        1: str= 'JAN'
        2: str= 'FEB'
        3: str= 'MAR'
        4: str= 'APR'
        5: str= 'MAY'
        6: str= 'JUN'
        7: str= 'JUL'
        8: str= 'AUG'
        9: str= 'SEP'
        10: str= 'OCT'
        11: str= 'NOV'
        12: str= 'DEC'
     endcase
     mnstr[i] = str
  endfor
  
end


;;decimal years to a new string
PRO yeartostr,year,str

  yr = fix(year)
  decimal = year-yr
  days = decimal*365
  month = 1
  while days gt 30 do begin
     case month of
        1: days-=31
        2: if yr mod 4 eq 0 then days-=29 else days-=28
        3: days-=31
        4: days-=30
        5: days-=31
        6: days-=30
        7: days-=31
        8: days-=31
        9: days-=30
        10: days-=31
        11: days-=30
        12: days-=31
     endcase
     month++
  endwhile
  print,yr,month,days

END



PRO wrapper

  print,'ARE YOU SURE?!??'
  stop
  readcol,'~/Dropbox/Orbit_MCMC/Orbit_test_sample.txt',name,format='a'
  name = strtrim(name,2)
  for i = 0,n_elements(name)-1 do begin
     formatex,name[i],/includebad,adderr=0
  endfor
  
END

;; just for objects in the binary sample
PRO wrapper2

  readcol,'~/Dropbox/Orbit_MCMC/finished.txt',name,format='a'
  name = strtrim(name,2)
  for i = 0,n_elements(name)-1 do begin
     formatex,name[i],/update,/erradd
  endfor

END


;; outputs formatted 'extra' data
PRO formatex,name,add=add,update=update,includebad=includebad,erradd=erradd

  if n_elements(erradd) eq 0 then erradd=1
  if n_elements(includebad) eq 0 then includebad=0
  if n_elements(update) eq 0 then update = 0
  if n_elements(add) eq 0 then add = 0
  restore,'~/Dropbox/delfosse/binaries.dat'
  
  l = where(binaries.name eq strtrim(name,2))
  if name eq 'GJ99' then l = where(binaries.name eq 'PM_I02287+3215')
  if name eq 'LHS1070AB' then l = where(binaries.name eq 'GJ2005AB')
  if name eq 'Gl795AC' then l = where(binaries.name eq 'PM_I20396+0458E')
  if name eq 'EGCha' then l = where(binaries.name eq 'EG_Cha')
  if l[0] eq -1 then begin
     querysimbad,name,ra,dec    ;,/cfa
     gcirc,2,ra,dec,binaries.ra,binaries.dec,dist
     l = where(dist eq min(dist) and dist lt 40)
     if l[0] eq -1 then begin
        print,'Failed to find target',min(dist)
        return
     endif
  endif
  if n_elements(l) gt 1 then begin
     print,'warning, two entries for the same star'
     l = l[0]
  endif
  b = binaries
  seps = *b[l].seps*1d0
  SEPS_ERR = *b[l].SEPS_ERR*1d0
  refs = strtrim(*b[l].ref,2)
  epochs = *b[l].epoch
  pas = *b[l].pas*1d0
  PAS_ERR = *b[l].PAS_ERR
  SEPS_FLAG = *b[l].SEPS_FLAG
  PAS_FLAG = *b[l].PAS_FLAG
  refs = strtrim(refs,2)
  baddies = ['PUEO_TMP','NACO_TMP','NIRI_TMP','TYC2000b','TYC2000','TYC2002','AST2001','AST2000','Los2010','Hrt1996a','Hry1998','WDK2015','Hen1997','Bnu1980b','Lei1994','Lei2001a','HIP1997a','AST1998','Beu2004','Gki2004','Sgr2000','TYC2000a']
  ;;Sgr2000 = PUEO
  if includebad eq 0 then begin
     readcol,'baddies.txt',bd,format='a'
     baddies = [baddies,bd]
     bd = strtrim(bd,2)
  endif
  
  ;; baddies = ['Ism1992','Inn1988','HIP1997a','Tok1982b','Tok1982a','Bnu1980b','Hry1998','WDK2015','TYC2000b','TYC2000','TYC2002','AST2001','AST2000','Beu2004','Bwl2015','PUEO_TMP','NACO_TMP','NIRI_TMP',$
  ;;            'Los2010',$  ;; this reference is questionable at best (french)
  ;;            'Car1996a',$ ;; too few points
  ;;            'Hrt1996a',$ ;; hard to categorize
  ;;            'Krv2016',$  ;; this one is GJ65 only, hard to calibrate
  ;;            'SJE2016',$  ;; uncalibratable
  ;;            'Dae2007',$  ;; uncalibratable
  ;;            'Dae2009',$  ;; uncalibatable
  ;;            ;;'Lei1994',$  ;; uncalibratable
  ;;            ;;'Lei2001a'$  ;; uncalibratable
  ;;            'Law2008',$  ;; uncalibratable
  ;;            'Mtg2006',$  ;; NACO/CFHT
  ;;            'RDR2015',$  ;; uncalibratable? (probably)
  ;;            'Los2010',$  ;; uncalibatable
  ;;            'SJE2014'$   ;; uncalibratable
  ;;          ] ;;,'McA1983','Gii2012','Cpb1994'
  match2,baddies,refs,suba,subb
  pas_flag = strtrim(pas_flag,2)
  ll = where(seps gt 0 and pas gt 0 and subb eq -1 and (PAS_FLAG eq '' or refs eq 'Masking' or pas_flag eq 'q' or pas_flag eq ':' or pas_flag eq 'V'))
  qq = where(seps gt 0 and (pas le 0 or PAS_FLAG ne '') and pas_flag ne 'q' and pas_flag ne ':')
  if qq[0] ne -1 then print,n_elements(qq),' epochs removed because of a missing P.A.'
  rr = where(subb ne -1 and strpos(refs,'_TMP') eq -1 and seps gt 0)
  if rr[0] ne -1 then begin
     print,n_elements(rr),' epochs removed because they are baddies:'
     forprint,epochs[rr],seps[rr],pas[rr],refs[rr]
     print,'---'
  endif
  if ll[0] eq -1 then begin
     print,'no useful data'
     return
  endif
  seps = seps[ll]*1d0
  pas = pas[ll]*1d0
  seps_err = seps_err[ll]*1d0
  pas_err = pas_err[ll]*1d0
  epochs = epochs[ll]*1d0
  pas_flag = pas_flag[ll]
  seps_flag = seps_flag[ll]
  refs = refs[ll]

  ;; references missing errors but for which we have derived from separate sources
  ;; ref_add = ['Hry1998','AST2016','Hor2017','Tok2015c','Tok2016a','Brg2010','USN1988b','WSI2004a','WSI2004b','AST2000','Mason2018'] ;; ,'Jnn2012' WSI2001b kinda sucks
  ;; sep_err_add = [2.0/1000d0,3.0/1000d0,1.7/1000d0,8.0/1000d0,8.0/1000d0,5.0/1d3,0.16,40/1000d0,40/1000d0,2/1000d0,6/1000d0]
  ;; pa_err_add = [0.5,0.5,0.5,1.5,1.5,1.0,1.5,1.5,1.5,0.0,1.75]
  ;; for jj = 0,n_elements(ref_add)-1 do begin
  ;;    gg = where(ref_add[jj] eq refs)
  ;;    if gg[0] ne -1 then begin
  ;;       seps_err[gg] = sep_err_add[jj]
  ;;       pas_err[gg] = pa_err_add[jj]
  ;;    endif
  ;; endfor
  
  ;; references where the measurement error needs to be increased
  ;; ref_add = ['Koh2012','Sef2008','Krv2016','Sgr2000','Jod2013','Jnn2014','Jnn2012','Bag2013','Bag2006b']
  ;; sep_err_add = [3.5/1000d0,0.0,0.0,7.0/1000d0,5./1000d0,5./1000d0,7./1000d0,4.0/1000d0,5.0/1000d0]
  ;; pa_err_add = [0.0,0.5,1.5,2.0,0.0,1.0,1.0,1.0,2.0]
  ;; for jj = 0,n_elements(ref_add)-1 do begin
  ;;    gg = where(ref_add[jj] eq refs)
  ;;    if gg[0] ne -1 then begin
  ;;       seps_err[gg] = sqrt(seps_err[gg]^2.+sep_err_add[jj]^2.)
  ;;       pas_err[gg] = sqrt(pas_err[gg]^2.+pa_err_add[jj]^2.)
  ;;    endif
  ;; endfor

  pa_sys_err = 0.0;0.5
  sep_sys_err = 0.0;2.5
  print,name
  qq = where(pas_err lt 0 or finite(pas_err) eq 0 or pas_err gt 100)
  if qq[0] ne -1 then pas_err[qq] = 0.0
  qq = where(seps_err lt 0 or finite(seps_err) eq 0 or seps_err gt 500)
  if qq[0] ne -1 then seps_err[qq] = 0.0

  ll = where(seps gt 0.0 and strpos(refs,'NIRC2') eq -1)
  if ll[0] eq -1 then begin
     print,'no data'
     return
  endif
  seps = seps[ll]
  pas = pas[ll]
  seps_err = seps_err[ll]
  pas_err = pas_err[ll]
  epochs = epochs[ll]
  pas_flag = pas_flag[ll]
  seps_flag = seps_flag[ll]
  refs = refs[ll]


  if name eq 'Gl795AC' then begin
     ll = where(seps lt 0.18)
     seps = seps[ll]
     pas = pas[ll]
     seps_err = seps_err[ll]
     pas_err = pas_err[ll]
     epochs = epochs[ll]
     pas_flag = pas_flag[ll]
     seps_flag = seps_flag[ll]
     refs = refs[ll]
  endif

  ;; check for uniqueness
  keep = intarr(n_elements(epochs))
  for iii = 0,n_elements(epochs)-1 do begin
     tmp = where(abs(epochs-epochs[iii]) lt 0.001 and abs(seps-seps[iii]) lt 0.001)
     if n_elements(tmp) gt 1 then begin
        if min(tmp) eq iii then keep[iii] = 1
     endif else keep[iii] = 1
     ;;if abs(epochs[iii]-2017.8246) lt 0.01 then stop
  endfor
  ll = wherE(keep eq 1)
  seps = seps[ll]
  pas = pas[ll]
  seps_err = seps_err[ll]
  pas_err = pas_err[ll]
  epochs = epochs[ll]
  pas_flag = pas_flag[ll]
  seps_flag = seps_flag[ll]
  refs = strtrim(refs[ll],2)


  
  ;;BIG ERROR, MASKING FILES ARE NOT HANDLED PROPERLY BECAUSE MISSING
  ;; ALSO WHAT DO WE DO FOR MISSING REFERENCES??
  ;; here is where we add the systematic error terms. These are stored
  ;; in a textfile
  if erradd eq 1 then begin
     readcol,'reference_catalog.txt',r,set,syssep,syspa,format='a,a,d,d',/silent
     syssep/=1000d0
     r = strtrim(r,2)
     for i = 0,n_elements(refs)-1 do begin
        if strtrim(refs[i],2) ne 'Masking' then begin
           l = where(strtrim(refs[i],2) eq strtrim(r,2))
           if l[0] eq -1 then begin
              ;; for now leave these as is??
              l = where(strtrim(set,2) eq strmid(refs[i],0,3))
              l = l[0]
              if l[0] ne -1 then print,'bad match? ',refs[i],set[l]
           endif
           if l[0] ne -1 then begin
              seps_err[i] = sqrt(seps_err[i]^2.+syssep[l]^2.)
              pas_err[i] = sqrt(pas_err[i]^2.+syspa[l]^2.)
           endif else begin
              print,'failed to find reference: ',refs[i]
              
           endelse
        endif
     endfor
     ll = wherE(seps_err lt 500 and pas_err lt 500)
     if ll[0] eq -1 then return 
     seps = seps[ll]
     pas = pas[ll]
     seps_err = seps_err[ll]
     pas_err = pas_err[ll]
     epochs = epochs[ll]
     pas_flag = pas_flag[ll]
     seps_flag = seps_flag[ll]
     refs = strtrim(refs[ll],2)
  endif

  if erradd eq 1 then begin
     qq = where(seps_err le 0 or pas_err le 0)
     if qq[0] ne -1 then begin
        print,'??'
        forprint,seps,pas,seps_err,pas_err,epochs,refs
        stop
     endif
  endif

     
  if ll[0] ne -1 then begin
     jday = JBEPOCH(epochs,/b,/to_day)
     DAYCNV, jday, YR, MN, DAY, HR
     round = where(hr gt 12.)
     if round[0] ne -1 then day[round]+=1
     getmnstr,mn,str
     daystr = string(day,format="(I02)")+'-'+str+'-'+string(yr,format="(I4)")
     mas = seps*1000d0
     e_mas = seps_err*1000d0
     e_pas = pas_err
     pas = pas
     refs = refs
     pflg = pas_flag
     flip = where(pflg eq 'q' and pas lt 180d0)
     if flip[0] ne -1 then pas[flip]+=180d0
     flip = where(pflg eq 'q' and pas gt 180d0)
     if flip[0] ne -1 then pas[flip]-=180d0
     ;;qq = where(e_mas lt 0.1) & if qq[0] ne -1 then e_mas[qq] = 15
     ;;qq = where(e_pas le 0.01) & if qq[0] ne -1 then e_pas[qq] = 2

     ;; old = where(strpos(daystr,'198') ne -1 or strpos(daystr,'197') ne -1)
     ;; if old[0] ne -1 then begin
     ;;    e_mas[old]*=3.0
     ;;    e_pas[old]*=3.0
     ;; endif
     ;; old = where(strpos(daystr,'199') ne -1)
     ;; if old[0] ne -1 then begin
     ;;    e_mas[old]*=1.5
     ;;    e_pas[old]*=1.5
     ;; endif
     ;; if name eq 'HD15285' then begin
     ;;    comment,daystr,'-DEC-1989'
     ;;    comment,daystr,'-APR-1991'
     ;;    comment,daystr,'-NOV-2000'
     ;;    comment,daystr,'-OCT-1978'
     ;; endif
     ;; if name eq 'GJ1005' then begin
     ;;    daystr[where(daystr eq '03-OCT-2002')] = '#03-OCT-2002'
     ;;    if name eq 'GJ1005' then daystr[where(refs eq 'Hry1998')] = '#'+daystr[where(refs eq 'Hry1998')] ;; duplicate with ast
     ;;    if name eq 'GJ1005' then daystr[where(daystr eq '30-NOV-2015')] = '#30-NOV-2015'                 ;;Tok2016a
     ;;    if name eq 'GJ1005' then daystr[where(daystr eq '22-JUN-2008')] = '#22-JUN-2008'                 ;;poor point
     ;; endif
     ;; if name eq 'GJ1245' then begin
     ;;    daystr[where(daystr eq '26-MAY-2005')] = '#26-MAY-2005'
     ;;    daystr[where(daystr eq '18-SEP-1986')] = '#18-SEP-1986'
     ;;    daystr[where(daystr eq '14-JUN-1987')] = '#14-JUN-1987'
     ;;    daystr[where(daystr eq '07-JUN-2008')] = '#07-JUN-2008'
     ;; endif
     ;; if name eq 'GJ4287' then begin
     ;;    daystr[where(daystr eq '02-APR-1991')] = '#02-APR-1991'
     ;; endif
     ;; if name eq 'GJ234' then begin
     ;;    daystr[where(daystr eq '32-DEC-1981')] = '#31-DEC-1981'
     ;;    daystr[where(daystr eq '05-OCT-1982')] = '#05-OCT-1982'
     ;;    daystr[where(refs eq 'Lmp2013')] = '#'+daystr[where(refs eq 'Lmp2013')]
     ;;    daystr[where(refs eq 'Bag2007b')] = '#'+daystr[where(refs eq 'Bag2007b')]
     ;;    daystr[where(refs eq 'Bag2006b')] = '#'+daystr[where(refs eq 'Bag2006b')]
     ;; endif
     ;; if name eq 'GJ1210' then begin
     ;;    pas[where(daystr eq '19-JUN-2008')]+=180d0                       ;; AB flip
     ;;    pas[where(daystr eq '24-JUN-2010')]+=180d0                       ;; AB flip
     ;;    pas[where(daystr eq '04-OCT-1982')]+=180d0                       ;; AB flip
     ;; endif
     ;; if name eq 'GJ3010' then begin
     ;;    daystr[where(daystr eq '07-AUG-2001')] = '#07-AUG-2001'
     ;;    daystr[where(daystr eq '08-JAN-2012')] = '#08-JAN-2012'
     ;; endif
     ;; if name eq 'Gl84' then daystr[where(daystr eq '03-JUL-2013')] = '#03-JUL-2013'
     ;; if name eq 'Gl84' then daystr[where(daystr eq '01-OCT-2002')] = '#01-OCT-2002'
     ;; if name eq 'HIP11542' then begin
     ;;    daystr[where(refs eq 'Doc2006i')] = '#'+daystr[where(refs eq 'Doc2006i')]                                  ;;'duplicate point'
     ;;    daystr[where(strpos(daystr,'AUG-1986') ne -1)] = '#'+daystr[where(strpos(daystr,'AUG-1986') ne -1)]        ;;junk
     ;;    daystr[where(refs eq 'Masking')] = '#'+daystr[where(refs eq 'Masking')]
     ;; endif
     ;; if name eq 'Gl125' then begin
     ;;    daystr[where(daystr eq '04-DEC-2003' and refs eq 'Bag2005')] = '#04-DEC-2003'                      ;;'duplicate point'
     ;;    daystr[where(daystr eq '04-OCT-2001' and refs eq 'Bag2005')] = '#04-OCT-2001'                      ;;'duplicate point'
     ;; endif
     ;; if name eq 'Gl804' then begin
     ;;    daystr[where(daystr eq '29-OCT-2004')] = '#'+daystr[where(daystr eq '29-OCT-2004')]
     ;;    daystr[where(daystr eq '12-JUN-2006')] = '#12-JUN-2006'
     ;; endif
     ;; if name eq 'Gl616.2' then begin
     ;;    comment,daystr,'-MAR-2007'
     ;;    comment,daystr,'-JUN-2007'
     ;;    pas[where(refs eq 'Bla1987')]+=180d0
     ;; endif
     ;; if name eq 'GJ340' then begin
     ;;    pas[where(daystr eq '26-APR-1999')]+=180d0                      ;; AB flip
     ;;    pas[where(daystr eq '20-DEC-2004')]+=180d0                      ;; AB flip
     ;;    daystr[where(daystr eq '07-JAN-1983')] = '#'+ daystr[where(daystr eq '07-JAN-1983')]
     ;;    daystr[where(daystr eq '24-JAN-1983')] = '#'+ daystr[where(daystr eq '24-JAN-1983')]
     ;;    daystr[where(daystr eq '20-JAN-1984')] = '#'+ daystr[where(daystr eq '20-JAN-1984')]
     ;;    ;;daystr[where(strpos(refs,'Hrt') ne -1)] = '#'+ daystr[where(strpos(refs,'Hrt') ne -1)]
     ;;    ;;daystr[where(strpos(refs,'WSI') ne -1)] = '#'+ daystr[where(strpos(refs,'WSI') ne -1)]
     ;;    ;;daystr[where(strpos(refs,'Pru') ne -1)] = '#'+ daystr[where(strpos(refs,'Pru') ne -1)]
     ;;    ;;daystr[where(strpos(refs,'Hrt2000a') ne -1)] = '#'+ daystr[where(strpos(refs,'Hrt2000a') ne -1)]
     ;;    ;;daystr[where(strpos(refs,'Pru') ne -1)] = '#'+ daystr[where(strpos(refs,'Pru') ne -1)]
     ;;    ;;e_pas+=50
     ;;    ;;e_mas+=100
     ;; endif
     ;; if name eq 'Gl54' then daystr[where(daystr eq '17-OCT-2008')] = '#17-OCT-2008'
     ;; if name eq 'Gl54' then daystr[where(daystr eq '08-OCT-2014')] = '#08-OCT-2014'
     ;; ;;if name eq 'GJ661' then daystr[where(refs eq 'WDK2015')] = '#'+daystr[where(refs eq 'WDK2015')]
     ;; if name eq 'GJ661' then daystr[where(refs eq 'Bag1989a')] = '#'+daystr[where(refs eq 'Bag1989a')]
     ;; if name eq 'Gl330' then pas[where(daystr eq '16-JAN-2014')]+=180d0 ;; AB flip
     ;; if name eq 'GJ65' then begin
     ;;    daystr[where(daystr eq '04-AUG-2010')]='#04-AUG-2010'                                         ;; AB flip
     ;;    daystr[where(daystr eq '28-NOV-1990')]='#28-NOV-1990'                                         ;; AB flip
     ;;    daystr[where(refs eq 'RDR2015')] = '#'+daystr[where(refs eq 'RDR2015')]
     ;; endif
     ;; if name eq 'GJ748' then daystr[where(refs eq 'AST1998')]='#'+daystr[where(refs eq 'AST1998')]    ;; duplicate
     ;; if name eq 'GJ748' then daystr[where(refs eq 'Mcy1983')]='#'+daystr[where(refs eq 'Mcy1983')]    ;; shit
     ;; if name eq 'GJ748' then daystr[where(refs eq 'AST2001')]='#'+daystr[where(refs eq 'AST2001')]    ;; shit
     ;; if name eq 'HD239960' then daystr[where(refs eq 'Los2010')]='#'+daystr[where(refs eq 'Los2010')] ;; shit
     ;; if name eq 'GJ2005' then begin
     ;;    daystr[where(refs eq 'Lei1994')]='#'+daystr[where(refs eq 'Lei1994')]   ;; junk
     ;;    daystr[where(refs eq 'Lei2001a')]='#'+daystr[where(refs eq 'Lei2001a')] ;; junk
     ;;    e_mas[where(refs eq 'Lei2001a')] = sqrt(e_mas[where(refs eq 'Lei2001a')]^2.+5.0^2.0)
     ;;    e_pas[where(refs eq 'Lei2001a')] = sqrt(e_pas[where(refs eq 'Lei2001a')]^2.+2.0^2.0)
     ;; endif
     ;; if name eq 'Gl310' then begin
     ;;    daystr[where(refs eq 'RDR2015')] = '#'+daystr[where(refs eq 'RDR2015')]
     ;; endif
     ;; if name eq 'GJ22' then daystr[where(refs eq 'TYC2002')]='#'+daystr[where(refs eq 'TYC2002')]   ;; junk
     ;; if name eq 'GJ22' then daystr[where(refs eq 'WSI2006a')]='#'+daystr[where(refs eq 'WSI2006a')] ;; junk
     ;; if name eq 'GJ22' then daystr[where(refs eq 'Jod2013')]='#'+daystr[where(refs eq 'Jod2013')]   ;; junk
     ;; if name eq 'Gl747' then begin
     ;;    daystr[where(refs eq 'Bag1991b')]='#'+daystr[where(refs eq 'Bag1991b')]                      ;; junk
     ;;    daystr[where(refs eq 'McA1987b')]='#'+daystr[where(refs eq 'McA1987b')]                      ;; junk
     ;;    pas[where(refs eq 'Bag1989a' and daystr ne '19-MAY-1987')]-=180d0
     ;;    pas[where(daystr eq '19-MAY-1987')]+=180d0
     ;;    pas[where(refs eq 'Bla1987')]-=180d0
     ;;    ;;e_pas[where(refs eq 'Sgr2000')] = sqrt(e_pas[where(refs eq 'Sgr2000')]^2.+0.5^2.0)
     ;;    ;;e_mas[where(refs eq 'Sgr2000')] = sqrt(e_mas[where(refs eq 'Sgr2000')]^2.+3.0^2.0)
     ;; endif
     if name eq 'HD283646' then pas[where(daystr eq '10-OCT-1998')]-=180d0
     if name eq 'GJ3412' then pas[where(daystr eq '27-OCT-2007')]-=180d0
     if name eq 'GJ3412' then daystr[where(daystr eq '12-NOV-2000')]='#'+daystr[where(daystr eq '12-NOV-2000')]
     if name eq 'GJ3421' then pas[where(daystr eq '20-JAN-2008')]+=180d0
     if name eq 'GJ3421' then  daystr[where(daystr eq '08-JAN-2012')]='#'+daystr[where(daystr eq '08-JAN-2012')]
     if name eq 'Gl469' then daystr[where(daystr eq '12-JUN-2009')]='#'+daystr[where(daystr eq '12-JUN-2009')]
     if name eq 'GJ600' then pas[where(daystr eq '21-JUN-2008' or daystr eq '01-MAY-2007' or daystr eq '12-FEB-2004')]-=180d0
     if name eq 'HD152751' then pas[where(daystr eq '15-JUL-2008' or daystr eq '12-JUN-2006' or daystr eq '01-JUL-2001'$
                                          or daystr eq '24-AUG-1999' or daystr eq '20-FEB-2000' or daystr eq '15-APR-2000'$
                                          or daystr eq '23-JUL-1998' or daystr eq '01-JUL-2015')]+=180d0
     if name eq 'HD152751' then  pas[where(daystr eq '04-APR-2002' or daystr eq '22-JUN-2000')]-=180
     if name eq 'HD152751' then e_pas[where(refs eq 'Jnn2012')] = sqrt(e_pas[where(refs eq 'Jnn2012')]^2.+0.5^2.0)
     if name eq 'Gl473' then  pas[where(daystr eq '11-MAR-2012' or daystr eq '14-JUN-1986')]+=180
     if name eq 'Gl473' then e_pas[where(refs eq 'Bla1987')] = sqrt(e_pas[where(refs eq 'Bla1987')]^2.+4.0^2.0)
     if name eq 'GJ3454' then begin
        daystr[where(refs eq 'Law2008')] = '#'+daystr[where(refs eq 'Law2008')]
        daystr[where(refs eq 'Hen1997')] = '#'+daystr[where(refs eq 'Hen1997')]
        ;; daystr[where(ref eq 'Jnn2014')] = '#'+daystr[where(ref eq 'Jnn2014')]
     endif
     if name eq 'Gl185' then begin
        daystr[where(daystr eq '19-SEP-2002')] = '#'+daystr[where(daystr eq '19-SEP-2002')] ;;redundant
     endif
     
     ;;if name eq 'GJ600' then  

     ;;ll = where(refs ne 'Masking')
     ;;if ll[0] ne -1 then begin
     ;;   e_pas[ll] = sqrt(e_pas[ll]^2.+pa_sys_err[ll]^2.)
     ;;   e_mas[ll] = sqrt(e_mas[ll]^2.+sep_sys_err[ll]^2.)
     ;;endif

     if add eq 0 and update eq 0 then begin
        forprint,daystr+string(9b)+string(mas,format="(D7.2)")+string(9b)+string(pas,format="(D6.2)")+string(9b)+string(e_mas,format="(D6.2)")+string(9b)+string(e_pas,format="(D6.2)")+string(9b)+refs+string(9b)+pflg,/nocomment,textout='extra_info/astrom-'+name+'AB-sp.txt',/silent
        print,n_elements(mas),'Epochs added'
     endif else begin
        if file_test('extra_info/astrom-'+name+'AB-sp.txt') then begin
           readcol,'extra_info/astrom-'+name+'AB-sp.txt',daystr_old,mas_old,pas_old,e_mas_old,e_pas_old,refs_old,format='a,d,d,d,d,a',/silent
           refs = strtrim(refs,2)
           refs_old = strtrim(refs_old,2)
           tmp1 = daystr_old
           tmp2 = daystr
           for jjj = 0,n_elements(tmp1)-1 do tmp1[jjj] = strreplace(tmp1[jjj],'#','')
           for jjj = 0,n_elements(tmp1)-1 do tmp1[jjj] = strreplace(tmp1[jjj],'#','')
           for jjj = 0,n_elements(tmp2)-1 do tmp2[jjj] = strreplace(tmp2[jjj],'#','')
           for jjj = 0,n_elements(tmp2)-1 do tmp2[jjj] = strreplace(tmp2[jjj],'#','')
           match2,strtrim(tmp2,2),strtrim(tmp1,2),suba,subb
           ll = where(suba eq -1)
           if ll[0] eq -1 then print,'No additional epochs added' else print,n_elements(ll),' additional epochs added'
           if add eq 1 then begin
              daystr = [daystr_old,daystr[where(suba eq -1)]]
              mas = [mas_old,mas[where(suba eq -1)]]
              pas = [pas_old,pas[where(suba eq -1)]]
              e_mas = [e_mas_old,e_mas[where(suba eq -1)]]
              e_pas = [e_pas_old,e_pas[where(suba eq -1)]]
              refs = [refs_old,refs[where(suba eq -1)]]
              
              
              if n_elements(pflg) lt n_elements(refs) then begin
                 pflg = [pflg,'']
              endif
           endif
           if update eq 1 then begin ;; erase all the old measurements and errors and replace them with better ones
              ;; preserve comments
              for iii = 0,n_elements(daystr)-1 do begin
                 ll = (where(refs[iii] eq refs_old and (daystr[iii] eq daystr_old or '#'+daystr[iii] eq daystr_old)))[0]
                 if ll[0] eq -1 then print,'Missing:',daystr[iii],refs[iii] else $
                    if strpos(daystr_old[ll],'#') ne -1 then daystr[iii] = '#'+daystr[iii]
              endfor
              
              ;; preserve P.A (flips)
              for iii = 0,n_elements(daystr)-1 do begin
                 ll = (where(refs[iii] eq refs_old and (daystr[iii] eq daystr_old or '#'+daystr[iii] eq daystr_old)))[0]
                 if ll[0] ne -1 then begin
                    pas[iii] = pas_old[ll]
                    mas[iii] = mas_old[ll]
                 endif
              endfor

           endif
        endif else begin
           print,'no existing file found',name
           return
        endelse
     endelse
     forprint,daystr+string(9b)+string(mas,format="(D7.2)")+string(9b)+string(pas,format="(D6.2)")+string(9b)+string(e_mas,format="(D6.2)")+string(9b)+string(e_pas,format="(D6.2)")+string(9b)+refs+string(9b)+pflg,/nocomment,textout='extra_info/astrom-'+name+'AB-sp.txt',/silent
  endif else begin
     print,'Failed to output any data'
  endelse
  
END
