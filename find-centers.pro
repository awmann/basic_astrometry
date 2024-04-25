
drive = 'Vali'

if n_elements(object) eq 0 then begin
   print,'Provide an object'
   return
endif   
;;object = '2MASSJ22240821+1728466';SSW J0920122+351742';'GJ416'


print
find_keck_data, object, filename=files, fail=fail
;;for i = 0,n_elements(files)-1 do files[i] = strreplace(files[i],'/Volumes/Loki/','/Volumes/Vali/')
for i = 0,n_elements(files)-1 do files[i] = strreplace(files[i],'/Volumes/Vali/','/Volumes/Loki/')
for i = 0,n_elements(files)-1 do files[i] = strreplace(files[i],'//','/')

;;AWM edit, checking against the blacklist before running this code
readcol, 'blacklist.txt', blackfiles, f='(a)', /silent, comment='#'
for i = 0,n_elements(blackfiles)-1 do blackfiles[i] = strreplace(blackfiles[i],'/Volumes/Vali/','/Volumes/Loki/')
for i = 0,n_elements(blackfiles)-1 do blackfiles[i] = strreplace(blackfiles[i],'//','/')
match2,files,blackfiles,suba,subb
help,where(suba ne -1)
print,'blacklist files: ',n_elements(where(suba ne -1))
files = files[where(suba eq -1)]
;; end AWM edit

if fail eq 1 then begin
   print,'no files found, stopping'
   stop
endif

nfiles = n_elements(files)
wclear
window, xs=750, ys=700

xguess = -1
yguess = -1

for i=0,nfiles-1 do begin

   fdecomp, files[i], disk__, dir__, name__, qual__
   outfile = 'centers/'+name__+'.xy.txt'
   
   tmpim = readfits(strreplace(files[i],'gfk',drive),/silent)
   ;if n_elements(tmpim[*,0]) lt 1024 or n_elements(tmpim[0,*]) lt 1024 then begin
   ;   old = tmpim
   ;   tmpim = reform_nirc2(tmpim)
   ;   tmpim[where(tmpim le 0)] = median(old)
   ;endif
   print
   print, files[i]
   spawn,'echo '+files[i]+' | pbcopy'
   count = 0
   happy = 0
   while happy ne 1 do begin

      ;; check for previous work
      if happy eq 0 and file_test(outfile) then begin
         ;; display image
         loadct,0,/silent
         display2, asinh_stretch(tmpim), /tv
         lincolr,/silent
         ; read in past data
         readcol, outfile, xprev, yprev, /silent
         oplot, xprev, yprev, psym=4, col=12, syms=3, thick=2
         happy = getyn('satisified with existing centers file?')
         if happy eq 0 then $
            spawn,'mv '+outfile+' '+outfile+'.bad'
         if happy eq 1 then begin
            xobj = xprev
            yobj = yprev
         endif
      endif

      ; check for previous clicks
      if happy eq 0 and xguess[0] ne -1 and yguess[0] ne -1 then begin
         happy=1
         
         ; display image
         loadct,0,/silent
         display2, asinh_stretch(tmpim), /tv
         lincolr,/silent
         ; centroid on previous guesses
         for k=0,total(finite(xguess))-1 do begin
            cntrd, tmpim, xguess[k], yguess[k], tmpxgc, tmpygc, 8.0
            if tmpxgc ne -1 and tmpygc ne -1 then begin
               xguess[k] = tmpxgc
               yguess[k] = tmpygc
            endif
         endfor
         oplot, xguess, yguess, psym=4, col=5, syms=3, thick=2
         happy = getyn('satisified with these positions?')
         if happy eq 1 then begin
            xobj = xguess
            yobj = yguess
         endif
      endif

      ; give up: just click everywhere
      if happy eq 0 then begin
         xobj = -1
         yobj = -1
         count = 0
         done_clicking = 0
         while done_clicking eq 0 do begin
            ; display image
            loadct,0,/silent
            display2, asinh_stretch(tmpim), /tv
            lincolr,/silent
            oplot, [xobj], [yobj], psym=4, col=3, syms=3, thick=2
            ; get new position
            cursor, tmpx, tmpy & wait,0.2
            ; if the user clicked on something then try to centroid there
            if tmpx gt 0 then begin
               print, count, tmpx, tmpy, f='(i4,2f10.3,$)'
               cntrd, tmpim, tmpx, tmpy, tmpxc, tmpyc, 6.0
               if tmpxc lt 0 or tmpyc lt 0 then begin
                  tmpxc = tmpx
                  tmpyc = tmpy
               endif else begin
                  ; and display the result, zoomed in
                  loadct,0,/silent
                  display2, asinh_stretch(tmpim[(tmpxc-20)>0:(tmpxc+20)<(n_elements(tmpim[*,0])-1),(tmpyc-20)>0:(tmpyc+20)<(n_elements(tmpim[0,*])-1)]), /tv
                  lincolr,/silent
                  oplot, [20], [20], psym=4, syms=5, thick=4, col=3
                  wait,0.5
                  loadct,0,/silent
                  display2, asinh_stretch(tmpim), /tv
                  lincolr,/silent
               endelse
               print, tmpxc, tmpyc, f='("   -->   ",2f10.3)'
               ; save this centroid
               xobj = [ tmpxc, xobj ]
               yobj = [ tmpyc, yobj ]
               count++
            endif else begin
               if tmpy gt 0 then begin
                  ; keep all clicks so far
                  nobj = n_elements(xobj)-1
                  xobj = xobj[0:nobj-1]
                  yobj = yobj[0:nobj-1]
                  done_clicking = 1
                  happy = 1
               endif else begin
                  ; remove last click
                  xobj = xobj[1:*]
                  yobj = yobj[1:*]
                  count--
               endelse    
            endelse
         endwhile
      endif
   endwhile
   
   ; determine the brightness order
   flux = xobj*0.
   for k=0,n_elements(xobj)-1 do begin
      flux[k] = total(tmpim[round(xobj[k])-1:round(xobj[k])+1,round(yobj[k])-1:round(yobj[k])+1],/nan)
   endfor
   
   ; output results
   forprint, text=outfile, xobj, yobj, flux/max(flux), files[i]+strarr(n_elements(xobj)), /silent, subset=sort(-flux), f='(2f10.3,x,f10.6,x,a)', /nocomment

   xguess = xobj
   yguess = yobj
   
   ;stop
      
endfor


end
