# Instructions for basic reduction

## Introduction

This guide will give you the simplest possible reduction, using the default behavior of the pipeline. It uses esorex, not the reflex GUI. So you do the reduction one step and a time, and you have the opportunity to inspect the data and understand what is going on. Also, at any stage you can change the parameters, or apply custom steps of your own. For example, you can improve the sky subtraction of each OB by median-subtracting the IFUs at each wavelengths (if you only care about line detections, and less about the actual fluxes), or re-calibrate the positional shifts of the dithering patterns and individual OBs by measuring the position of helper stars, and apply those shifts to the final reconstructed cubes. This is not covered here, there is too much already ;)

Now are you prepared to loose a full day or more? Then let's go.

NB: You can use these scripts freely. In exchange, a citation to this webpage or Schreiber et al. (in prep.) is always welcome!

## A. Prepare your computer

1) Download the ESO pipeline:
ftp://ftp.eso.org/pub/dfs/pipelines/kmos/kmos-kit-1.3.14.tar.gz
Extract it into some temporary directory and run the "./install_pipeline" script.

2) Make sure that the "bin" directory of the kmos pipeline is in your PATH so that the executables can be reached from anywhere.

3) Define the Bash variable KMOS_CALIB_DIR to point to the directory of the KMOS pipeline that contains the static calibration data. This directory must contain files like "kmos_oh_spec_hk.fits". For me it is (put in .bashrc or .profile or whatever you use in Mac):
```bash
export KMOS_CALIB_DIR="/home/cschreib/programming/kmos/calib/kmos-1.3.14/cal/"
```

4) Install the QFitsView tools:
http://www.mpe.mpg.de/~ott/dpuser/qfitsview.html

5) Some of the scripts I introduce below are written in C++ and use the phy++ library I developed during my PhD. To install it, follow the following instructions:
- Download this: https://github.com/cschreib/phypp/archive/master.tar.gz
- Extract it into some temporary directory.
- Make a new directory there called "build" and go there.
- Make sure the development files for cfitsio are installed on your computer. If you have troubles with this step, tell me.
- Call `cmake ../` (NB: you need to install CMake for this to work)
- Then call `make`, and `sudo make install`.
- Make sure to follow the last instructions that were printed in the terminal ("sourcing" the file ".phypprc")

6) Attached is a tar.gz archive containing the scripts I use. Extract it somewhere, wherever. We will copy the scripts one by one in their correct directories later on.

## B. Prepare the data

1) Uncompress the data given by ESO.
```bash
for f in *.Z; do gzip -d $f; done
```

2) With a terminal go into this directory, and copy the the attached phy++ program "rename.cpp". You can run it with the command below; it will extract the purpose of each FITS file in the directory and rename it accordingly.
```bash
ophy++ rename ./
```

  Names:
  - `*sci.fits`: science images (for both A and B)
  - `*acq.fits`: acquisition images (bright stars used to position the telescope before the science runs)
  - `*sky.fits`: sky images (never used)
  - `*object-sky-std-flux.fits`: flux calibrators images
  - `*dark.fits`, `*flat-off.fits`, `*flat-lamp.fits`, `*wave-off.fits`, `*wave-lamp.fits`, `*flat-sky.fits`: various calibration images
  (you can ignore the files that start with M.KMOS.*)

3) Now manually group the files into new directories based 1) on the time at which they were taken, and 2) based on the category:
- `*sci.fits`, `*sky.fits`, `*acq.fits`, `*NL.txt` and `*.xml` files go into `sci-XX`
- `*object-sky-std-flux.fits` go into `calib-std-XX`
- the rest go into `calib-XX`

"XX" is a number starting from "01" and increasing for each group. Basically, you want to group together the calibration files that were taken on the same day, and science images that belong to a single OB. I decide this based on the name of the files, which contains the date of observation. Be careful that, since observations are done in the night, sometime the date can change in the middle of an OB or calibration set ;) Similarly, sometimes two OBs can be executed on the same day. Be sure to check the hours and minutes, and that each OB contains the right number of "sci" images (e.g., for us we have an AABAA pattern, so 5 "sci" images per OB). Each OB is typically preceded by several "acq" frames, sometimes "sky" as well. Standard stars calibration sets are always made of 4 "std-flux" images.

Other calibration sets can contain various types of images, not always the same. Usually this is 5 "dark", 3 "flat-off", 18 "flat-lamp", 1 "wave-off" and 6 "wave-lamp". Often there are 10 "dark", but the pipeline can only handle 5, so I only keep the 5 latest and place the earliest 5 into a sub-directory "unused" (this is important else you can get errors in the pipeline). In addition, sometimes you have 4 "flat-sky" a bit before the calibration set. I don't know how standard this is, so I did not automatize this part.

For one of our KMOS run, if I do this, I get 6 OBs, 4 calibration sets and 5 standard star calibration sets.

## C. Reduce the calibration data

1) Then the boring part... Reducing the calibration data. Make sure that all the files that you manipulated in the previous step are stored in a separate and safe directory. Say, `/home/cschreib/data/kmos/`. To keep things clean, create a new directory, for example `/home/cschreib/data/kmos/reduced/`. I'll call this the "working directory", and the rest of the work will be done here in order to avoid touching the raw data by accident.

2) To make the task a bit simpler I have created a script to automatize most of the process. It is called "reduce.cpp". Copy it inside the working folder and compile it:
```bash
cphy++ optimize reduce.cpp
```

3) Now copy the attached "make_calib.sh" file into the working directory. Open it to make sure that :
- all the calibration sets are listed in the CALIBS variable
- the raw data directory is correct in RAW_DIR
- the value of GRATING matches your observations (the format is `XXX`, where "X" is the name of the band in which you observe; by default is it the K band, so `GRATING=KKK`; in our some other of our own KMOS programs we used H+K, so `GRATING=HKHKHK`; you get the point).

Then run it.
```bash
chmod +x make_calib.sh
./make_calib.sh
```

For each calibration set, this script will create several ".sof" files needed by the pipeline. These files list the raw data files (plus some other stuff) that the pipeline will reduce. It also creates a "reduce.sh" script inside each directory, that can be called to actually do the reduction (see next step). These "reduce.sh" scripts will appear over and over in the rest of the guide: normally you don't have to modify them, but there is no magic, all the calls to the pipeline are there. You can tweak the parameters there if you wish.

NB: You will notice that "make_calib.sh" will print some warnings. Typically:
```bash
note: no 'flat-sky' frame in /home/cschreib/data/kmos/calib-01/
warning: no illumination correction in this set
```

Indeed, the people at ESO do the "flat-sky" calibration very rarely, which is needed for illumination correction (compensate for the fact that some IFUs receive a bit more light than others). This is an optional step in the science reduction however. In any case, this is normal, and I will show you later how to deal with this.

However if any other file is missing, the reduction cannot proceed and "make_calib.sh" will tell you so.

4) To actually reduce the data, you can go inside each "calib-XX" directory and run the "reduce.sh" script. Alternatively, if you want to reduce everything in one command, you can copy the attached "reduce_calib.sh" file into the working directory, run it and go get a coffee.
```bash
chmod +x reduce_calib.sh
./reduce_calib.sh
```

5) Once the calibration is finished, you can take a look at the calibration end products (open with DS9) to make sure that everything went fine. Below I give a short description of the main calibration files, and tell you what to expect based on my own experience. You don't have to look at this right now, there is a lot! But if something goes wrong or just if you want to understand what is going on, you should probably read this.

- badpixel_dark.fits: flag "hot" pixels in the detector with 0, and good pixels with 1. The file has 3 extensions, one for each detector. This "badpixel" map is only used internally by the calibration reduction, science frames use the badpixel maps presented below. Most of the pixels should be equal to 1, and many small spots should be flagged with 0.

- badpixel_flat_XXX.fits: same as above, but with 3x6 extensions. For each 3 detectors, this image contains the badpixel mask at 6 different telescope rotator angle (see the keyword "ESO PRO ROT NAANGLE"). Your science data uses a single rotator angle, and the pipeline will pick the closest badpixel mask in this file. In addition to "badpixel_dark.fits", these images also flag out the pixels that are not illuminated by any IFU. The resulting file should look like "badpixel_dark.fits", but with 14x8 vertical stripes.

- master_dark.fits: "zero flux" image of each detector, when the telescope is not observing anything (not even the sky). The FITS extensions of these files are the same as for "badpixel_dark.fits" maps (3 extensions). The maps should be fairly homogeneous, with values close to 0.

- master_flat_XXX.fits: "uniform illumination" image of each detector with the IFUs, when the telescope is observing a calibration lamp. The FITS extensions of these files are the same as for "badpixel_flat_XXX.fits" (3x6 extensions). The maps should be fairly homogeneous, with values close to 1.

- xcal_XXX.fits, ycal_XXX.fits and lcal_XXX.fits: These are probably the most important calibration files. The format is the same as "badpixel_flat_XXX.fits" (3x6 extensions). Combined, these three files tell you the sky and wavelength coordinate of each pixel of each detector. That is: to which IFU the pixel belongs, then within this IFU, to which spatial pixel (NB: there are 14x14 spatial pixels in an IFU), and finally what wavelength of the spectrum. The IFU identifier is given as the decimal part of "xcal" and "ycal". For example, if a pixel value is "1400.5", then it means that the pixel belongs to the IFU 5. IFUs are numbered from 1 to 8 within each detector. In these two files, the integral part (i.e., "1400" in the example above) is the sky position in milliarcseconds, relative to the center of the IFU. The format of "lcal" is straightforward, the value of each pixel is simply the wavelength in microns.

Sometimes, you may see that there is a huge vertical "gap" in one of the calibration images, for example in "badpixel_flat", "master_flat", "xcal", "ycal" and "lcal". This simply means that one of the IFU was disabled when they did the calibration (and probably, also for your observation, but it may be good to check).

If something goes wrong, e.g. if one of the files above is missing or weird, open the "reduce.sh" script and launch each reduction step one at a time. Then look at the "esorex.log" file to see if it is telling you something. It happened to us in a previous KMOS program that one of the calibration set was bugged, and could not be reduced. In this case the solution is to forget about this calibration set and use another one for the observations that depended on it.

## D. Reduce the standard stars for absolute flux calibration

1) Now that the calibration is fully reduced, we can reduce the standard stars to get absolute flux calibration. You may not care of the actual flux for your science case, but I think it can help reduce the RMS (not sure though). This is short anyway, and it will prepare you for the real science reduction later on since the procedure is very similar. So let's go! Copy the attached "make_stdstar.sh" script into the working directory and open it with your text editor.

Now, as for "make_calib.sh", make sure that the values of `RAW_DIR` and `GRATING` are correct. Then you see two things. The first is a Bash function that I use as a shortcut ("reduce_wrapper"). Don't touch it. Then there is a list of lines like:
```bash
reduce_wrapper 01 calib=[../calib-01,../calib-03] # 2015-09-27 08:28
```

Here is how it works. The line must start with "reduce_wrapper". Then you must give the ID of the standard star calibration set, i.e., "01" here. If you did what I told you in point (B.3) above, this is the "01" in "calib-std-01". You will want one line for each "calib-std-XX" in your raw data. Finally, the most tricky part is `calib=[...]` (the comment after "#" is optional, just for you: I usually write the observation date here)

Within the "[...]" you want to list the calibration sets that the pipeline will use to reduce the spectrum of the standard stars. Normally, you would only give one calibration set, e.g., "calib-01". But as you may have seen in point (C.3), not all calibration sets have the illumination correction. The trick is that we can use the illumination correction of *another* calibration set. So the idea is to put in the "[...]" first the "main" calibration set, the one that is closest in time (either before or after) to the observations you want to reduce, and second the other calibration set from which the illumination correction will be taken. It is important that the "main" set is given first in the list. If the "main" calibration set already contains the illumination calibration, there is no need to give a second set. Got it? I'll give you an example to make it clear.

Example: calib-std-01, observed from 2015-12-26T06:50:53.933 to 2015-12-26T06:54:58.577.
The closest calibration set is calib-01, observed 2015-12-26T09:34:15.697 to 2015-12-26T10:36:37.706, i.e., right after.
So we write:
```bash
reduce_wrapper 01 calib=[../calib-01]
```
(nb: don't forget the "../" in front of the calibration set)

However, "calib-01" doesn't include illumination correction. If you noted down the result of step (C.3), you will find (for example) that in fact only "calib-03" has it. So we add it to the list:
```bash
reduce_wrapper 01 calib=[../calib-01,../calib-03]
```

And that's it! Now do the same for all the other standard star sets.

2) Once you have figured out the associations of all the calibration stars, you can run the script to create the pipeline ".sof" files and reduction scripts:
```bash
chmod +x make_stdstar.sh
./make_stdstar.sh
```

3) ... and reduce this. As for the calibration, either you go into each "calib-std-XX" directory and run the "reduce.sh" script, or you copy the attached "reduce_stdstar.sh" in the working directory and run it.
```bash
chmod +x reduce_stdstar.sh
./reduce_stdstar.sh
```

4) If you wish you can inspect the result of this reduction, to see if everything went well. In each "calib-std-XX" directory the pipeline produces several files, including the full KMOS spectral cube, the extracted spectrum, and the continuum image. Use QFitsView to inspect those.

## E. A short tutorial to inspect a cube with QFitsView

1) Launch QFitsView, go to `File -> Open` and navigate to the directory where the cube is located.

2) In the file list, select the file containing the data cube, but do not open it yet. First untick the box "All Extensions" (if it was ticked) and select "1" in the "Image Extension" list (the default, "primary", is empty). Then open it.

3) This will show the cube of the 1st IFU in the file. Unfortunately, for the standard stars, only 3 IFUs are actually used, so you may have to open a different extension to actually see something.

4) Now QFitsView shows you one image in the middle of the cube, the spectral element 1024 (this number can be found in a small edit box in the tool bar at the top of the window). What you usually want to do first is to display the continuum image, i.e., the sum of the cube on the wavelength axis. You can do that by choosing "Average" or "Median" in the drop down list to the right of "1024" (by default it is set to "Single"). For the standard star it doesn't change much, since the S/N is very high even for each spectral element, but for your galaxies you usually only detect them in the continuum.

5) Then you can move your mouse over the image. At the bottom of the window you should see the spectrum change in real time: this is the spectrum of the pixel currently below your mouse. Since sources are usually spread over several pixels (because they are extended or because of the seeing), you usually want to average multiple neighboring pixels together. You can do that by changing the "Single pixel" (close to the spectrum) into "Circular", and increase the radius "R1" to 2 or 3 pixels. This is simply aperture photometry :) You will see that now the shape of your aperture displayed under the mouse pointer, and the spectrum is now averaged over all these pixels, enhancing the S/N.

6) If you want to save the spectrum to analyze it with another tool, like IDL, put your mouse where you want to extract the spectrum, then right-click and choose "Save spectrum as..." (ASCII or FITS, as you prefer).

7) You can also "lock" the spectrum, i.e., freeze the position of the aperture, by right-clicking and choosing "Lock position"

8) If you have identified a line and want to see if it is real or not, a good thing to do is to look at the "line image" of your galaxy. To do that you want to change "Average" (or "Median") into "Linemap" in the top toolbar. The current wavelength that is displayed is shown in the spectrum with a gray vertical bar. You can change the wavelength either by specifying the ID of the spectral element (default is 1024, and it goes from 1 to 2048), or by moving the slider just below.

9) Once you are at the right wavelength, maybe the S/N is not very high because you are displaying only a single spectral element, while the line is spectrally extended (it has a width). You can increase the S/N by averaging multiple spectral elements at once. To do that, go to "Options -> Cube display...".  Choose "Line map", and in "Channels", increase the second value of "Central frame". The default is 1, so you only display one spectral element. You may want to try 3 or 6 or 10. This should give a better image of the line spatial profile.

## F. Reduce individual OBs

1) Now that you have been through the standard star calibration, this will seem easy :) Copy the "make_sci.sh" script into the working directory and open it with your text editor. The format is the same as for "make_stdstar.sh", and the only difference is that you also have to provide which standard star calibration set to use for the flux calibration:
```bash
reduce_wrapper sci-01 calib=[../calib-01,../calib-03] stdstar=../stdstar-01 # 2015-12-26T06:04:36.909
```

Again, choose the standard star that is observed the closest to your science data. It is preferable if the standard star is reduced with the same calibration sets (`calib=[...]`).

2) Make the script executable and run it.
```bash
chmod +x make_sci.sh
./make_sci.sh
```

3) As usual, you can go into each "sci-XX" directory and run the "reduce.sh" script, or copy the "reduce_sci.sh" script into the working directory and run it.
```bash
chmod +x reduce_sci.sh
./reduce_sci.sh
```

4) Now each OB has been reduced individually, with the sky subtracted. The pipeline can recognize automatically if you have put science targets in both the A and B configurations, so you don't have anything special to do. In a next step we will merge these OBs into a single thing, but for now it is important to make sure that the reduction was successful (no bug or crash). Just go into each "sci-XX" directory: there should be `N` FITS files, corresponding to each of the `N` dither positions you specified in your OBs. Each FITS file contains 48 extensions, two per IFU (one for the flux and the other for the uncertainty). Disabled IFUs are also included here, but the extension does not contain any data. You don't need to check them one by one (see next point), just make sure the files are there.

## G. Check for the detection of helper targets

1) Next you want to make sure that your "helper" targets can be seen in the continuum image for each OB. These are stars of magnitude 21-19, as recommended in the manual. Additionally, some times we chose to observe a z=0 galaxy with a Pashen-alpha line: it will show both a continuum and line detection in each OB, so it can be used to check the wavelength calibration and how the line is affected by the reconstruction.

NB: At this stage there are some other things we can do to improve the reduction, like improving the sky subtraction and astrometry, but we will see that some other time. For now we will just check that the helper targets are well detected. To do these checks, you could use QFitsView and open each reduced cube one by one. This is tedious... Instead you can follow the procedure below, which I find more convenient.

2) First you are going to need the names of your helper targets, as given in the KARMA file when you prepared the OBs. To find out the names of all your targets, go into, e.g., "sci-01" and run the following command on any of the FITS images:
```bash
fitshdr sci_reconstructed_KMOS.2015-12-26T05\:56\:49.884-sci.fits | grep -E "ARM[0-9]+ NAME"
```

3) This will print the list of the IFUs "ARMx" (where "x" is the IFU ID) and the corresponding targets. From there, identify the name of your helper targets. You may have to repeat the above command for another OB, if you have multiple target lists.

4) Now copy the attached phy++ program "extract_ifu.cpp" into the working directory and compile it:
```bash
cphy++ optimize extract_ifu.cpp
```

This tool extracts the IFUs given their target name, since I couldn't find a way to do it with the pipeline.

5) Then copy the attached "collapse_helpers.sh" script into the working directory, and open it with your text editor. Make sure that: a) `GRATING` has the correct value, b) `HELPERS` contains the full list of the helper target names you want to reduce, with format "[name1,name2,name3,etc]" (no spaces!) and c) all your science OBs are listed into `SCIS`.

6) Make it executable and run it:
```bash
chmod +x collapse_helpers.sh
./collapse_helpers.sh
```

7) Then for each of the "sci-XX" directories, run the following command:
```bash
ds9 sci-XX/helpers/*.fits
```

This will open DS9 and display the continuum images of all your helper targets in the OB, for each of the 4 dither positions. Offsets are normal and expected, for now just make sure that they are detected. Non-detections can be caused by a number of factors, including wrong acquisition, wrong calibration, etc. Weak detections can be caused by bad seeing.

## H. Combine all OBs into master cubes

1) Copy the attached "make_combine.sh" script into the working directory. Then make it executable and run it.
```bash
chmod+x make_combine.sh
./make_combine.sh
```

2) This creates a new directory called "sci-master". Go there, and run the "reduce.sh" script.
The pipeline now combines all OBs and the dither patterns, taking into account the dithering position offset. It will create two files per target, one is the actual data cube, and the other is the exposure map which tells you how many DITs were observed at each spatial position of the IFU (it is not that useful, since it is not taking into account flagged pixels etc).

3) Now you can play with QFitsView to inspect your final data :)

4) At this stage, I like to take a look at the continuum images of all the targets, to have an idea of what is going on. To automatically generate these images, create a new directory inside "sci-master", for example "continuum", go there and copy the attached "make_collapsed.sh" script. Open it with your text editor and make sure that the grating is correct. Then make it executable and run it, then run "reduce.sh". For each data cube, it will create an "*_img_cont.fits" file with the continuum image. You can then open them all at once in DS9:
```bash
ds9 *_img_cont.fits
```
