# Keck astrometry fitting code

This assumes basic reduced data from NIRC2 Keck. It will output astrometry. Note that this code/ does not correct for differences between NIRC2 and Gaia or other zero-point offsets. 

I have assumed you are on my computer (Wuotan), and hence the paths follow that convention. update accordingly. I usually work in a folder called "BINFIT". 

A basic run looks like this:

    .r update-weather-savfile.pro
	This code assumes the cfht weather files are in BINFIT
    find_keck_data,/update_cat
	This will record new Keck data /Volumes/Baldur/NIRC2_reduced/
	This will take a while to run
	Note that this will only run on Wuotan
    .r build_binary
    list_finished

That was just book keeping. The core of the code you probably want is the rest of it:

find-centers.pro
     This will have the user find the centers by eye. 

binfit-keck.pro
     Run this right after find-centers. You will need to adjust sf_fit =0,1 as appropriate. That is to say, if you want to use starfinder or not. Starfinder usually is superior, but it sometimes fails catastrophically. 
     You can take the output from binfit-keck as is, but it's just in X-Y space. 

compute-sep-pa.pro
     this inherits the obj name (assuming you ran this right after binfit) and converts X/Y to separation and PA, using the distortion solution for NIRC2. 