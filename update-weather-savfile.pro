;; get the data from here:
;; http://mkwc.ifa.hawaii.edu/archive/wx/cfht/
;; then simply fill out the information below

file_old = 'cfht-wx.2020.dat'
file_new = 'cfht-wx.2021.dat'

spawn,'rm mko.tmp mko.tmp3'
spawn,'cp '+file_old+' mko.tmp'
spawn,'sed ''s/$/  NaN   NaN/'' mko.tmp > mko.tmp3'
readcol, 'mko.tmp3', yr_old, mo_old, dy_old, hr_old, mi_old, wspd_old, wdir_old, temp_old, relh_old, pres_old, $
         f='(d,d,d,d,d,d,d,d,d,d)'
jdcnv, yr_old, mo_old, dy_old, 10d + hr_old+mi_old/60d, jd_old
mjd_mko_old = jd_old-2400000.5d

spawn,'rm mko.tmp mko.tmp3'
spawn,'cp '+file_new+' mko.tmp'
spawn,'sed ''s/$/  NaN   NaN/'' mko.tmp > mko.tmp3'
readcol, 'mko.tmp3', yr_new, mo_new, dy_new, hr_new, mi_new, wspd_new, wdir_new, temp_new, relh_new, pres_new, $
         f='(d,d,d,d,d,d,d,d,d,d)'
jdcnv, yr_new, mo_new, dy_new, 10d + hr_new+mi_new/60d, jd_new
mjd_mko_new = jd_new-2400000.5d


restore, 'mko-weather.sav';, mjd_mko, wspd_kts, wdir_deg, temp_c, relh_pct, pres_mb

w = where( mjd_mko gt max(mjd_mko_old) or mjd_mko eq mjd_mko_old or mjd_mko eq mjd_mko_new, n )
if n gt 0 then $
   remove, w, mjd_mko, wspd_kts, wdir_deg, temp_c, relh_pct, pres_mb

mjd_mko  = [ mjd_mko, mjd_mko_old, mjd_mko_new ]
wspd_kts = [ wspd_kts , wspd_old, wspd_new  ]
wdir_deg = [ wdir_deg , wdir_old, wdir_new  ]
temp_c   = [ temp_c   , temp_old, temp_new  ]
relh_pct = [ relh_pct , relh_old, relh_new  ]
pres_mb  = [ pres_mb  , pres_old, pres_new  ]

save, file='mko-weather.sav', mjd_mko, wspd_kts, wdir_deg, temp_c, relh_pct, pres_mb


end
