;+
; FIND_KECK_DATA
;
;
;
;
;
;
;
;
;
;-
pro find_keck_data, ra_in, de_in, $          ;; ra/dec of object to search for
                    update_cat=update_cat, $ ;; ??
                    show_files=show_files, $ ;;
                    search_rad=search_rad, $
                    filename=filename, $
                    name_in=name_in, $
                    drive=drive, $ ;; e.g., Vali/Kraus_nirc2/Reduced
                    masking_check=masking_check, $
                    silent=silent, $
                    fail = fail

  if n_elements(silent) eq 0 then silent = 0
  fail = 0
  if not keyword_set(search_rad) then search_rad=60. ; arscec
  if not keyword_set(name_in) then name_in=''
  ;;if not keyword_set(drive) then drive='Loki'
  if not keyword_set(drive) then drive='Baldur'

  base = 'ldif'

  imglist = '~/Dropbox/BINFIT/'+base+'-nirc2-img.txt'
  catfile = '~/Dropbox/BINFIT/'+base+'-nirc2-cat.txt'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CREATE/UPDATE MASSIVE CATALOG OF IMAGES
  if keyword_set(update_cat) then begin

     print,'updating '+base+' catalog...'

     if file_test(catfile) then begin
        readcol, catfile, imfile_, datestr_, ra_, de_, objname_, filter_, f='(a,a,f,f,a,a)',/silent
        imfile_ = '/Volumes/'+drive+strmid(imfile_,strpos(imfile_[0],'/NIRC2_reduced/'))
     endif else begin
        imfile_ = ''
        datestr_ = ''
        ra_ = !values.f_nan
        de_ = !values.f_nan
        objname_ = ''
        filter_ = ''
     endelse

     print,'Generating imglist, '+imglist
     ;;spawn,'find /Volumes/'+drive+'/NIRC2_reduced/?????/????/ -type f -name ''N2.????????.?????.LDIF.fits'' > '+imglist
     spawn, 'wc '+imglist

     readcol, imglist, imfile, f='(a)',/silent
     print,imfile[-1]

     nimg = n_elements(imfile)
     datestr = strarr(nimg)+'0000-00-00'
     objname = strarr(nimg)+'ERROR_NO_NAME'
     filter = strarr(nimg)+'ERROR_NO_NAME'
     ra = strarr(nimg)+'NaN'
     de = strarr(nimg)+'NaN'
     count = 0l
     for i=0,nimg-1 do begin
        w = where( imfile_ eq imfile[i], n )

        if n eq 0 then begin
           hd = headfits(imfile[i])
           datestr[i] = strc(sxpar(hd,'DATE-OBS'))
           objname[i] = strc(sxpar(hd,'OBJECT'))
           filter[i] = strc(sxpar(hd,'FILTER'))
           ra[i] = strc(sxpar(hd,'RA'))
           de[i] = strc(sxpar(hd,'DEC'))
           if ra[i] eq '' or ra[i] eq '####Error###' or ra[i] eq '0' then ra[i] = 'NaN'
           if de[i] eq '' or de[i] eq '####Error###' or de[i] eq '0' then de[i] = 'NaN'
           if datestr[i] eq '' or datestr[i] eq '####Error###' or datestr[i] eq '0' then datestr[i] = '0000-00-00'
           if objname[i] eq '' or objname[i] eq '####Error###' or objname[i] eq '0' then objname[i] = 'ERROR_NO_NAME'
           if filter[i] eq '' or filter[i] eq '####Error###' or filter[i] eq '0' then filter[i] = 'ERROR_NO_NAME'
           count++
        endif else begin
           datestr[i] = datestr_[w[0]]
           ra[i] = ra_[w[0]]
           de[i] = de_[w[0]]
           objname[i] = objname_[w[0]]
           filter[i] = filter_[w[0]]
        endelse
     endfor
     print, strc(count)+' new files added'
     w = where( strpos(ra,'null') ne -1 or strpos(de,'null') ne -1, n )
     if n gt 0 then begin
        ra[w] = 'NaN'
        de[w] = 'NaN'
     endif
     w = where( float(ra) lt 0. or float(ra) gt 360. or float(de) lt -90. or float(de) gt 90., n )
     if n gt 0 then begin
        ra[w] = 'NaN'
        de[w] = 'NaN'
     endif
     forprint, text=catfile, imfile, datestr, ra, de, objname, filter, f='(a-'+strc(max(strlen(imfile)))+',3x,a-'+strc(max(strlen(datestr)))+',3x,a-'+strc(max(strlen(ra)))+',3x,a-'+strc(max(strlen(de)))+',3x,a-'+strc(max(strlen(objname)))+',3x,a-'+strc(max(strlen(filter)))+')', subset=multisort(datestr,imfile)
     print, 'testing (valid lines should be same number above)...'
     readcol, catfile, test1, test2, test3, test4, test5, test6, f='(a,a,f,f,a,a)',/silent
     ra = float(ra)
     de = float(de)
     save, file=strreplace(catfile,'.txt','.sav'), catfile, imfile, datestr, ra, de, objname, filter

  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;readcol, catfile, imfile, datestr, ra, de, objname, f='(a,a,f,f,a)', /silent
  restore, strreplace(catfile,'.txt','.sav')

  for i=0,n_elements(ra_in)-1 do begin
     if type(ra_in[i]) eq 7 then begin
        if ra_in[i] eq 'GJ3193B' then begin
           querysimbad, 'GJ3193', tmpr, tmpd, found=tmp ;;,/cfa
        endif else begin
           if ra_in[i] eq 'LHS1070AB' then begin
              querysimbad, 'GJ2005A', tmpr, tmpd, found=tmp ;;,/cfa
           endif else begin
              if ra_in[i] eq '2MASSJ04595855-0333123' or  ra_in[i] eq 'PMJ04599-0333E' then begin
                 tmpr = 74.994011499999999
                 tmpd = -3.5533923000000001
                 tmp = 1
              endif else begin
                 if ra_in[i] eq 'Gl795AC' then begin
                    tmpr = 309.907318
                    tmpd = 4.971921
                    tmp = 1
                 endif else begin
                    querysimbad, ra_in[i], tmpr, tmpd, found=tmp ;;,/cfa
                 endelse
              endelse
           endelse
        endelse
        if tmp eq 0 then $
           message, 'input name not found in SIMBAD: '+ra_in[i]
        if silent eq 0 then print, 'INPUT TARGET: '+string(tmpr,f='(f09.5)')+' '+string(tmpd,f='(f+09.5)')+' '+ra_in[i]
     endif else begin
        tmpr = ra_in[i]
        tmpd = de_in[i]
        if name_in[0] ne '' and silent eq 0 then $
           print, 'INPUT TARGET: '+string(tmpr,f='(f09.5)')+' '+string(tmpd,f='(f+09.5)')+' '+name_in[i] $
        else $
           if silent eq 0 then print, 'INPUT TARGET: '+string(tmpr,f='(f09.5)')+' '+string(tmpd,f='(f+09.5)')
     endelse

     gcirc, 1, ra/15., de, tmpr/15., tmpd, dist
     w = where( dist lt search_rad, n )
     if n eq 0 then begin
        if silent eq 0 then print,'--- no data ---'
        fail = 1
     endif else begin
        maskstr = strarr(n)
        if keyword_set(masking_check) then begin
           for k=0,n-1 do begin
              if strpos(filter[w[k]],'hole') ne -1 then begin
                 maskstr[k] = 'masking'
              endif
           endfor
        endif
        if keyword_set(show_files) then begin
           w2 = multisort(datestr[w],imfile[w])
           for k=0,n-1 do $
              if silent eq 0 then  print, k, imfile[w[w2[k]]], datestr[w[w2[k]]], ra[w[w2[k]]], de[w[w2[k]]], objname[w[w2[k]]], maskstr[w2[k]], f='(i4,x,a-'+strc(max(strlen(imfile[w])))+',x,a-'+strc(max(strlen(datestr[w])))+',x,a-'+strc(max(strlen(ra[w])))+',x,a-'+strc(max(strlen(de[w])))+',x,a-'+strc(max(strlen(objname[w])))+',x,a)'
        endif else begin
           udate = datestr[w[uniq(datestr[w],sort(datestr[w]))]]
           for j=0,n_elements(udate)-1 do begin
              w2 = where( datestr[w] eq udate[j], n2 )
              maskstr2 = ''
              if keyword_set(masking_check) then $
                 maskstr2 = '('+string(total(maskstr[w2] ne ''),f='(i3)')+' masking)   '
              if silent eq 0 then print, udate[j], n2, maskstr2, imfile[w[w2[0]]], f='(a-10,3x,i3," files   ",a,a-'+strc(max(strlen(imfile[w])))+')'
           endfor
        endelse
        filename = imfile[w]
     endelse
     print
  endfor


end


;;FUNCTION strc,string
;;  return,strtrim(string,2)
;;END
