# Keck astrometry fitting code

This assumes basic reduced data from NIRC2 Keck. It will output astrometry. Note that this code/ does not correct for differences between NIRC2 and Gaia or other zero-point offsets. The core of the code is StarFinder (http://www.bo.astro.it/StarFinder/paper6.htm) or a triple-Gaussian. (sf_fit=1 or 0). You can also fit with an empirical template, which yields similar results. To do this, you need to replace the model wrapped in mpfit to a template then add in the PSF location, stretch, and scale factors as free parameters. 

I have assumed you are on my computer (Wuotan), and hence the paths follow that convention. update accordingly. I usually work in a folder called "BINFIT". 

A basic run looks like this:

	.r update-weather-savfile.pro
This code assumes the cfht weather files are in BINFIT

	find_keck_data,/update_cat
This will record new Keck data /Volumes/Baldur/NIRC2_reduced/

	.r build_binary
	list_finished
Lists the new data

That was just bookkeeping. The core of the code you probably want is the rest of it:

	find-centers.pro
This will have the user find the centers by eye. 

	binfit-keck.pro
Run this right after find-centers. You will need to adjust sf_fit =0,1 as appropriate. That is to say, if you want to use starfinder or not. Starfinder usually is superior, but it sometimes fails catastrophically. 
     You can take the output from binfit-keck as is, but it's just in X-Y space. 

	compute-sep-pa.pro
this inherits the obj name (assuming you ran this right after binfit) and converts X/Y to separation and PA, using the distortion solution for NIRC2. 


If you use this code, please cite the starfinder code:
https://arxiv.org/abs/astro-ph/0004101
as well as one of the following:
https://ui.adsabs.harvard.edu/abs/2019ApJ...871...63M/abstract
https://ui.adsabs.harvard.edu/abs/2016ApJ...817...80D/abstract
