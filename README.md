# Instructions for basic and advanced reduction

## Introduction

This guide will teach you how to obtain the simplest possible reduction of your VLT-KMOS data, using the default behavior of the pipeline. It uses esorex, not the reflex GUI. So you do the reduction one step and a time, and you have the opportunity to inspect the data and understand what is going on. Also, at any stage you can change the parameters, or apply custom steps of your own: the last sections describe some improvements over the standard reduction.

Now are you prepared to loose a full day or more? Then let's go. The whole guide may appear very long, but there are actually very few things for you to do; although probably much more to learn if you are new to near-IR spectropy.

NB: You can use these scripts freely. In exchange, a citation to this webpage or Schreiber et al. (in prep.) is always welcome!

## A. Prepare your computer

1) Download the ESO pipeline (be sure to pick up the latest version, the link below may be outdated):
ftp://ftp.eso.org/pub/dfs/pipelines/kmos/kmos-kit-1.3.17.tar.gz
Extract it into some temporary directory and run the "./install_pipeline" script.

2) Make sure that the "bin" directory of the kmos pipeline is in your PATH so that the executables can be reached from anywhere.

3) Define the Bash variable KMOS_CALIB_DIR to point to the directory of the KMOS pipeline that contains the static calibration data. This directory must contain files like "kmos_oh_spec_hk.fits". For me it is (put in .bashrc or .profile or whatever you use in Mac):
```bash
export KMOS_CALIB_DIR="/home/cschreib/programming/kmos/calib/kmos-1.3.17/cal/"
```

4) Install the QFitsView tools:
http://www.mpe.mpg.de/~ott/dpuser/qfitsview.html

5) Some of the scripts in this suite are written in C++ and use the phy++ library I developed during my PhD. To install it, follow the following instructions:
- Download this: https://github.com/cschreib/phypp/archive/master.tar.gz
- Extract it into some temporary directory. This should create a folder called `phypp-master`.
- Step inside this directory and create a new directory there, called "build". Go there with a terminal.
- Make sure the development files for cfitsio are installed on your computer. If you have troubles with this step, tell me.
- From within the `build` directory, call `cmake ../` (NB: you need to install CMake for this to work)
- Then call `make`, and `sudo make install`.
- Make sure to follow the last instructions that were printed in the terminal ("sourcing" the file ".phypprc")

## B. Prepare the data

1) Uncompress the data given by ESO.
```bash
for f in *.Z; do gzip -d $f; done
```

2) With a terminal go into this directory, and copy the phy++ program "rename.cpp". You can run it with the command below; it will extract the purpose of each FITS file in the directory and rename it accordingly.
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

Other calibration sets can contain various types of images, not always the same. Usually this is 5 "dark", 3 "flat-off", 18 "flat-lamp", 1 "wave-off" and 6 "wave-lamp". Often there are 10 "dark", but the pipeline can only handle 5, so I only keep the 5 latest and place the earliest 5 into a sub-directory "unused" (this is important else you can get errors in the pipeline). In addition, sometimes you have 4 "flat-sky" a bit before the calibration set.

I don't know how standard all of this is, so I did not automatize this part. Here is an exerpt of the directory structure for one of our own programs:
```bash
 + calib-01/
 - - KMOS.2014-10-05T20:56:06.900-dark.fits
 - - KMOS.2014-10-05T20:57:14.819-dark.fits
 - - KMOS.2014-10-05T20:58:22.630-dark.fits
 - - KMOS.2014-10-05T20:59:31.072-dark.fits
 - - KMOS.2014-10-05T21:00:38.882-dark.fits
 - - KMOS.2014-10-05T22:36:18.174-flat-sky.fits
 - - KMOS.2014-10-05T22:36:50.848-flat-sky.fits
 - - KMOS.2014-10-05T22:37:07.382-flat-sky.fits
 - - KMOS.2014-10-05T22:37:22.459-flat-sky.fits
 - - KMOS.2014-10-06T10:49:52.240-flat-off.fits
 - - KMOS.2014-10-06T10:50:02.071-flat-off.fits
 - - KMOS.2014-10-06T10:50:13.116-flat-off.fits
 - - KMOS.2014-10-06T10:53:24.259-flat-lamp.fits
 - - KMOS.2014-10-06T10:53:35.191-flat-lamp.fits
 - - KMOS.2014-10-06T10:53:46.146-flat-lamp.fits
 - - KMOS.2014-10-06T10:54:47.996-flat-lamp.fits
 - - KMOS.2014-10-06T10:54:58.914-flat-lamp.fits
 - - KMOS.2014-10-06T10:55:08.748-flat-lamp.fits
 - - KMOS.2014-10-06T10:56:10.491-flat-lamp.fits
 - - KMOS.2014-10-06T10:56:21.556-flat-lamp.fits
 - - KMOS.2014-10-06T10:56:32.503-flat-lamp.fits
 - - KMOS.2014-10-06T10:57:34.377-flat-lamp.fits
 - - KMOS.2014-10-06T10:57:45.316-flat-lamp.fits
 - - KMOS.2014-10-06T10:57:56.322-flat-lamp.fits
 - - KMOS.2014-10-06T10:58:58.075-flat-lamp.fits
 - - KMOS.2014-10-06T10:59:09.118-flat-lamp.fits
 - - KMOS.2014-10-06T10:59:20.035-flat-lamp.fits
 - - KMOS.2014-10-06T11:00:21.786-flat-lamp.fits
 - - KMOS.2014-10-06T11:00:32.700-flat-lamp.fits
 - - KMOS.2014-10-06T11:00:43.676-flat-lamp.fits
 - - KMOS.2014-10-06T11:46:21.039-wave-off.fits
 - - KMOS.2014-10-06T11:49:42.675-wave-lamp.fits
 - - KMOS.2014-10-06T11:50:46.487-wave-lamp.fits
 - - KMOS.2014-10-06T11:51:50.448-wave-lamp.fits
 - - KMOS.2014-10-06T11:52:54.252-wave-lamp.fits
 - - KMOS.2014-10-06T11:53:58.104-wave-lamp.fits
 - - KMOS.2014-10-06T11:55:01.959-wave-lamp.fits
 + calib-02/
 - - KMOS.2014-12-11T23:23:28.670-flat-sky.fits
 - - KMOS.2014-12-11T23:24:40.262-flat-sky.fits
 - - KMOS.2014-12-11T23:24:58.400-flat-sky.fits
 - - KMOS.2014-12-11T23:25:15.404-flat-sky.fits
 - - KMOS.2014-12-12T09:42:15.623-dark.fits
 - - KMOS.2014-12-12T09:43:23.428-dark.fits
 - - KMOS.2014-12-12T09:44:31.209-dark.fits
 - - KMOS.2014-12-12T09:45:39.011-dark.fits
 - - KMOS.2014-12-12T09:46:46.841-dark.fits
 - - KMOS.2014-12-12T15:36:46.045-flat-off.fits
 - - KMOS.2014-12-12T15:36:57.034-flat-off.fits
 - - KMOS.2014-12-12T15:37:07.968-flat-off.fits
 - - KMOS.2014-12-12T15:39:26.899-flat-lamp.fits
 - - KMOS.2014-12-12T15:39:36.698-flat-lamp.fits
 - - KMOS.2014-12-12T15:39:47.746-flat-lamp.fits
 - - KMOS.2014-12-12T15:40:49.458-flat-lamp.fits
 - - KMOS.2014-12-12T15:41:00.535-flat-lamp.fits
 - - KMOS.2014-12-12T15:41:11.458-flat-lamp.fits
 - - KMOS.2014-12-12T15:42:13.333-flat-lamp.fits
 - - KMOS.2014-12-12T15:42:24.247-flat-lamp.fits
 - - KMOS.2014-12-12T15:42:35.186-flat-lamp.fits
 - - KMOS.2014-12-12T15:43:37.078-flat-lamp.fits
 - - KMOS.2014-12-12T15:43:47.999-flat-lamp.fits
 - - KMOS.2014-12-12T15:43:58.966-flat-lamp.fits
 - - KMOS.2014-12-12T15:45:00.677-flat-lamp.fits
 - - KMOS.2014-12-12T15:45:11.711-flat-lamp.fits
 - - KMOS.2014-12-12T15:45:22.903-flat-lamp.fits
 - - KMOS.2014-12-12T15:46:24.647-flat-lamp.fits
 - - KMOS.2014-12-12T15:46:35.608-flat-lamp.fits
 - - KMOS.2014-12-12T15:46:45.410-flat-lamp.fits
 - - KMOS.2014-12-12T16:08:15.344-wave-off.fits
 - - KMOS.2014-12-12T16:11:36.775-wave-lamp.fits
 - - KMOS.2014-12-12T16:12:40.642-wave-lamp.fits
 - - KMOS.2014-12-12T16:13:44.439-wave-lamp.fits
 - - KMOS.2014-12-12T16:14:48.236-wave-lamp.fits
 - - KMOS.2014-12-12T16:15:52.025-wave-lamp.fits
 - - KMOS.2014-12-12T16:16:55.805-wave-lamp.fits
 + calib-std-01/
 - - KMOS.2014-10-10T09:18:24.913-object-sky-std-flux.fits
 - - KMOS.2014-10-10T09:19:54.063-object-sky-std-flux.fits
 - - KMOS.2014-10-10T09:21:19.090-object-sky-std-flux.fits
 - - KMOS.2014-10-10T09:22:44.392-object-sky-std-flux.fits
 + calib-std-02/
 - - KMOS.2014-11-25T04:06:57.978-object-sky-std-flux.fits
 - - KMOS.2014-11-25T04:08:25.027-object-sky-std-flux.fits
 - - KMOS.2014-11-25T04:09:50.356-object-sky-std-flux.fits
 - - KMOS.2014-11-25T04:11:16.333-object-sky-std-flux.fits
 + sci-01/
 - - KMOS.2014-11-25T05:35:39.767-acq.fits
 - - KMOS.2014-11-25T05:35:39.767.NL.txt
 - - KMOS.2014-11-25T05:35:39.767.xml
 - - KMOS.2014-11-25T05:38:19.499-acq.fits
 - - KMOS.2014-11-25T05:38:19.499.NL.txt
 - - KMOS.2014-11-25T05:38:19.499.xml
 - - KMOS.2014-11-25T05:43:42.076-sci.fits
 - - KMOS.2014-11-25T05:43:42.076.NL.txt
 - - KMOS.2014-11-25T05:43:42.076.xml
 - - KMOS.2014-11-25T05:54:19.911-sci.fits
 - - KMOS.2014-11-25T06:04:38.076-sci.fits
 - - KMOS.2014-11-25T06:14:56.200-sci.fits
 - - KMOS.2014-11-25T06:25:09.193-sci.fits
 + sci-02/
 - - KMOS.2014-11-25T07:09:51.615-acq.fits
 - - KMOS.2014-11-25T07:09:51.615.NL.txt
 - - KMOS.2014-11-25T07:09:51.615.xml
 - - KMOS.2014-11-25T07:11:28.902-acq.fits
 - - KMOS.2014-11-25T07:11:28.902.NL.txt
 - - KMOS.2014-11-25T07:11:28.902.xml
 - - KMOS.2014-11-25T07:16:38.243-sci.fits
 - - KMOS.2014-11-25T07:16:38.243.NL.txt
 - - KMOS.2014-11-25T07:16:38.243.xml
 - - KMOS.2014-11-25T07:27:17.464-sci.fits
 - - KMOS.2014-11-25T07:37:35.540-sci.fits
 - - KMOS.2014-11-25T07:47:53.040-sci.fits
 - - KMOS.2014-11-25T07:58:06.107-sci.fits
```

For this particular KMOS program (10 hours), I get in the end 10 OBs, 9 calibration sets and 13 standard star calibration sets (I show only two of each above, since the whole list would otherwise be very long). You may get more or less similar numbers, depending on how long was your program.

## C. Reduce the calibration data

1) Then the boring part... Reducing the calibration data. Make sure that all the files that you manipulated in the previous step are stored in a separate and safe directory. Say, `/home/cschreib/data/kmos/`. To keep things clean, create a new directory, for example `/home/cschreib/data/kmos/reduced/`. I'll call this the "working directory", and the rest of the work will be done here in order to avoid touching the raw data by accident.

2) To make the task a bit simpler I have created a script to automatize most of the process. It is called "reduce.cpp". Copy this file inside the working folder and compile it:
```bash
cphy++ optimize reduce.cpp
```

3) Now copy the "make_calib.sh" file into the working directory. Open it to make sure that :
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

4) To actually reduce the data, you can go inside each "calib-XX" directory and run the "reduce.sh" script. Alternatively, if you want to reduce everything in one command, you can copy the "reduce_calib.sh" file into the working directory, run it and go get a coffee.
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

1) Now that the calibration is fully reduced, we can reduce the standard stars to get absolute flux calibration. You may not care of the actual flux for your science case, but I think it can help reduce the RMS (not sure though). This is short anyway, and it will prepare you for the real science reduction later on since the procedure is very similar. So let's go! Copy the "make_stdstar.sh" script into the working directory and open it with your text editor.

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

3) ... and reduce this. As for the calibration, either you go into each "calib-std-XX" directory and run the "reduce.sh" script, or you copy the "reduce_stdstar.sh" in the working directory and run it.
```bash
chmod +x reduce_stdstar.sh
./reduce_stdstar.sh
```

4) I advise you to inspect the result of this reduction, to see if everything went well. In each "calib-std-XX" directory the pipeline produces several files, including the full KMOS spectral cube, the extracted spectrum, and the continuum image. Use QFitsView to inspect those (see next section for a short tutorial). In particular, take look at each of the 3 extensions of `std_image_XXX.fits`. These are collapsed images of the standard stars; in principle you should see a bright source at the center of the IFU.

In one of the calibration set we received, two out of three of these stars were almost out of the field of view of the IFU. The reduction went without telling us, but the flux calibration was completely off, increasing the noise in the final reduced science data by up to 20%. In some other cases the star could not even be seen in the IFU. The solution in all these cases was to ignore these faulty standard stars and use another set that was observed a couple of hours after.

## E. A short tutorial to inspect a cube with QFitsView

1) Launch QFitsView, go to `File -> Open` and navigate to the directory where the cube is located.

2) In the file list, select the file containing the data cube, but do not open it yet. First untick the box "All Extensions" (if it was ticked) and select "1" in the "Image Extension" list (the default, "primary", is empty). Then open it.

3) This will show the cube of the 1st IFU in the file. Unfortunately, for the standard stars, only 3 IFUs are actually used, so you may have to open a different extension to actually see something.

4) Now QFitsView shows you one image in the middle of the cube, the spectral element 1024 (this number can be found in a small edit box in the tool bar at the top of the window). What you usually want to do first is to display the continuum image, i.e., the sum of the cube on the wavelength axis. You can do that by choosing "Average" or "Median" in the drop down list to the right of "1024" (by default it is set to "Single"). For the standard star it doesn't change much, since the S/N is very high even for each spectral element, but for your galaxies you usually only detect them in the continuum.

5) Then you can move your mouse over the image. At the bottom of the window you should see the spectrum change in real time: this is the spectrum of the pixel currently below your mouse. Since sources are usually spread over several pixels (because they are extended or because of the seeing), you usually want to average multiple neighboring pixels together. You can do that by changing the "Single pixel" (close to the spectrum) into "Circular", and increase the radius "R1" to 2 or 3 pixels. This is simply aperture photometry :) You will see now the shape of your aperture displayed under the mouse pointer, and the spectrum is now averaged over all these pixels, enhancing the S/N.

6) If you want to save the spectrum to analyze it with another tool, like IDL, put your mouse where you want to extract the spectrum, then right-click and choose "Save spectrum as..." (ASCII or FITS, as you prefer).

7) You can also "lock" the spectrum, i.e., freeze the position of the aperture, by right-clicking and choosing "Lock position"

8) If you have identified a line and want to see if it is real or not, a good thing to do is to look at the "line image" of your galaxy. To do that you want to change "Average" (or "Median") into "Linemap" in the top toolbar. The current wavelength that is displayed is shown in the spectrum with a gray vertical bar. You can change the wavelength either by specifying the ID of the spectral element (default is 1024, and it goes from 1 to 2048), or by moving the slider just below.

9) Once you are at the right wavelength, maybe the S/N is not very high because you are displaying only a single spectral element, while the line is spectrally extended (it has a width). You can increase the S/N by averaging multiple spectral elements at once. To do that, go to "Options -> Cube display...".  Choose "Line map", and in "Channels", increase the second value of "Central frame". The default is 1, so you only display one spectral element. You may want to try 3 or 6 or 10. This should give a better image of the line spatial profile, and may also allow you to detect the continuum emission if you choose a large enough spectral window.

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

4) Now each OB has been reduced individually, with the sky subtracted. The pipeline can recognize automatically if you have put science targets in both the A and B configurations, so you don't have anything special to do. In a next step we will merge these OBs into a single thing, but for now it is important to make sure that the reduction was successful (no bug or crash). Just go into each "sci-XX" directory: there should be `N` FITS files, corresponding to each of the `N` exposures you specified in your OBs. Each FITS file contains 48 extensions, two per IFU (one for the flux and the other for the uncertainty). Disabled IFUs are also included here, but the extension does not contain any data. You don't need to check them one by one (see next point), for now just make sure the files are there.

If there is nothing, then the reduction failed. You can inspect the content of the `esorex_sci.log` file to figure out what happened and, if possible, fix the issue.

## G. Check for the detection of helper targets

1) Next you want to make sure that your "helper" targets can be seen in the continuum image for each OB. These are stars of magnitude H=21 (for DIT=600s), as informally recommended in the manual. Some times we chose instead (or in addition) to observe a z=0 galaxy with a Pashen-alpha line: it will show both a continuum and line detection in each OB, so it can be used to check the wavelength calibration and how the line is affected by the reconstruction. If you did not observe specific targets for this purpose, you can also use one of your science targets, provided it at least as bright in the NIR continuum and not too extended.

NB: At this stage there are some other things we can do to improve the reduction, like improving the sky subtraction and astrometry, but we will see that later on. For now we will just check that the helper targets are well detected. To do these checks, you could use QFitsView and open each reduced cube one by one. This is tedious... Instead you can follow the procedure below, which I find more convenient.

2) First you are going to need the names of your helper targets, as given in the KARMA file when you prepared the OBs. To find out the names of all your targets, go into, e.g., "sci-01" and run the following command on any of the FITS images:
```bash
fitshdr sci_reconstructed_KMOS.2015-12-26T05\:56\:49.884-sci.fits | grep -E "ARM[0-9]+ NAME"
```

3) This will print the list of the IFUs "ARMx" (where "x" is the IFU ID) and the corresponding targets. From there, identify the name of your helper targets. You may have to repeat the above command for another OB, if you have multiple target lists.

4) Now copy the C++ program "extract_ifu.cpp" into the working directory and compile it:
```bash
cphy++ optimize extract_ifu.cpp
```

This tool extracts the IFUs given their target name, since I couldn't find a way to do it with the pipeline.

5) Then copy the "collapse_helpers.sh" script into the working directory, and open it with your text editor. Make sure that: a) `GRATING` has the correct value, b) `HELPERS` contains the full list of the helper target names you want to reduce, with format "[name1,name2,name3,etc]" (no spaces!) and c) all your science OBs are listed into `SCIS`.

6) Make it executable and run it:
```bash
chmod +x collapse_helpers.sh
./collapse_helpers.sh
```

7) Then for each of the "sci-XX" directories, run the following command:
```bash
ds9 sci-XX/helpers/sci_reconstructed_*.fits
```

This will open DS9 and display the continuum images of all your helper targets in the OB, for each of the exposures. If you used dithering, offsets are normal and expected (try to match the WCS astrometry of the images), for now just make sure that they are detected. Non-detections can be caused by a number of factors, including wrong acquisition. We will see in section (J) how to deal with these cases (do not bother doing it right now). Weak detections can be caused by bad seeing.

9) You can also look at the combined images, merging together all the exposures of this particular OB:
```bash
ds9 sci-XX/helpers/combine_sci_reconstructed_*.fits

10) Using these combined images, another good test to do is to pick one helper target in particular and look at its continuum images in all the OBs of your program at onces. Assuming you chose the target whose name is `xxx`, run the following command from the working directory:
```bash
ds9 $(find | cat | sort | grep "helpers/combine_sci_reconstructed_xxx_img_cont.fits")
```

There, issues that you may notice concern the flux calibration and OH line subtraction. If you find that some of your OBs have much higher/lower noise and flux than the others, this may be indicative that the standard stars used for the flux calibration were improperly reduced. In this case you will want to check the photometric zero point computed by the pipeline and see if it deviates substantially from the other exposures.

You can also observe important background level variations, with some exposures having significantly negative or positive background. This is typical when OH line subtraction was imperfect, and I give a simplistic way to fix that later in section (I). But before you fix this, it is good to first perform a naive reduction without trying to fix anything, and this is what we do in the next section.

## H. Combine all OBs into master cubes

1) Copy the "fill_nan.cpp" file into the working directory and compile it:
```bash
cphy++ optimize fill_nan.cpp
```
This tool is used to flag out the wavelengths for which there is no observation. By default the pipeline puts a value of zero for these wavelengths, but these zeroes are included in the computation of the continuum images... So this tool fills these regions with NaN pixels, which are properly discarded by the pipeline.

2) Copy the "make_combine.sh" script into the working directory. Then make it executable and run it.
```bash
chmod+x make_combine.sh
./make_combine.sh
```

3) This creates a new directory called "sci-master" with the usual reduction script and SOF file. All the exposures of your program will be combined. If you have identified problematic exposures (for example because the helper targets could not be detected or strongly off-centered), you can remove the corresponding files from the "combine.sof" file. Then run the "reduce.sh" script.

The pipeline now combines all OBs and the dither patterns, taking into account the dithering position offset. It will create two files per target, one is the actual data cube, and the other is the exposure map which tells you how many exposures contribute to each spatial position of the IFU (it is not that useful, since it is not taking into account flagged pixels etc).

4) Now you can play with QFitsView to inspect your final data :)

5) At this stage, I like to take a look at the continuum images of all the targets, to have an idea of what is going on. To automatically generate these images, create a new directory inside "sci-master", for example "continuum", go there and copy the "make_collapsed.sh" script. Open it with your text editor and make sure that the grating is correct. Then make it executable and run it, then run "reduce.sh". For each data cube, it will create an "*_img_cont.fits" file with the continuum image. You can then open them all at once in DS9:
```bash
ds9 *_img_cont.fits
```

Using these images you can check further the sky position accuracy of the IFUs. However, depending on how bright your targets are, they may not appear at all in the continuum. Do not let that discourage you!

# I. Improve OH line subtraction

In some of our programs, the OH line subtraction is not optimal because our sky exposures are taken too far apart in time compared to the science exposures. The net result is that the sky lines are shifted in wavelength and/or intensity, and leave strong residuals. A simple empirical workaround for this issue is, for each IFU and each wavelength slice, to subtract the median of the pixel values. This is assuming the target is small compared to the IFU (for galaxies which are too large, part or most of their flux will be lost, so keep that in mind).

1) Copy the `median_sub.cpp` file into the working directory and compile it:
```bash
cphy++ optimize median_sub.cpp
```

2) Copy the `make_msub.sh` script into the working directory, make it executable and run it:
```bash
chmod +x make_msub.sh
./make_msub.sh
```

3) This create a copy of each reduced OB, and performs the median subtraction. The original OBs are always preserved, in case you are not satisfied with the result. Then you can come back to section H to re-reduce the OBs. To do so, just open the `make_combine.sh` script and modify the last line to:
```bash
reduce_wrapper sci-msub
```

Then run it and see if this improves the final data quality. If not, you should probably come back to the original OBs since the median subtraction can introduce biases in the flux measurements.

# J. Fix astrometry of individual exposures

The point of observing helper targets (see section G) is to make sure that the IFUs are properly centered on the targets in each exposure before collapsing them into a final cube. Most of the time it will be the case. However there is always the possibility that the acquisition step (done by the observer on site) was imperfect, either because of technical issues, or because the acquisition stars you provided are faulty (inaccurate position, or proper motion). This translates into positional shifts of the helper targets into the IFU. Below is the procedure to follow to correct for these shifts. Steps (1) to (10) should be done for each OB of your program.

1) For each OB, go into the `helpers` directory (which was created in section G). Choose one of the helper target (I'll call it `xxx`) and open the file `shifts_applied_xxx.txt` with your text editor. This file tells you which positional shifts the pipeline expects between all the exposures of this OB, knowing which dithering pattern you chose (if any). The first exposures always has zero shifts since it is used as a reference: the values of `x` and `y` are the centering difference between that exposure and the first one, given in pixels.

2) Copy the values of these shifts (omitting the first line with zero shifts) into a new file called `shifts.txt`. The format should be:
```bash
dx1 dy1
dx2 dy2
... ...
```
The file therefore contains two columns, with no header or comments, and one less lines than the number of exposures in the OB.

3) Open the `reduce.sh` script with your text editor, and locate the group of lines that start with the comment:
```bash
# Full xxx
```
Then modify the first line below to look like this:
```bash
esorex kmos_combine -method='user' -filename='shifts.txt' -cmethod='median' -name='xxx' combine.sof
```

4) This particular step is just a check to make sure you are doing everything correctly. You don't need to do it, and if you do, just do it for the first OB. Make a copy of the `combine_sci_reconstructed_xxx_img_cont.fits` file associated to your chosen helper target. Then run `reduce.sh` and make sure that the new `combine_sci_reconstructed_xxx_img_cont.fits` file shows the exact same image as the old one. Indeed, we applied manually the same shifts than the one calculated by the pipeline. If the images differ, then you have incorrectly written your `shifts.txt` file.

5) Now we will determine the real position shifts from the reduced images. Open the exposures with DS9:
```bash
ds9 *_img_cont-xxx.fits
```

6) Select the first frame, zoom to see it better, and change the contrast so that the helper target is clearly visible. Match these parameters to the other frames using  `Frame -> Match -> Frame -> WCS`, `Frame -> Match -> Scale` and `Frame -> Match -> Color bar`.

8) Now, in the first frame, draw a circular region centered on the helper target and copy this region into the other frames. If the exposures are correctly aligned with one another, the circular region should fall exactly on top of your helper target in all exposures (within one pixel or less). If this is the case, there is nothing to do and you can go to the next OB.

9) Else, you need to find the correct shifts. This is rather simple: first go to the first frame, double-click the circular region and note down its physical center coordinates (i.e., change `fk5` into `Physical`), I'll call them `x0` and `y0`. Then, for each of the other exposures, move the circular region to match the actual centroid of your helper target, and look at the resulting physical coordinates `xi` and `yi`. The shift is obtained simply by computing `xi - x0` and `y0 - yi` (yes, this is a convention of the pipeline). Use these values to replace the corresponding line in the `shift.txt` file.

10) Now do the procedure of step (4) above and compare the image of your helper target after and before applying the custom shifts. If you did well, the peak flux of the helper target should have increased, since its profile is less blurred by the position uncertainty. If this is not the case, then you may have computed the wrong shifts, or used the wrong order for the exposures, or maybe the shift corrections were too small to have any measurable impact.

11) Once you have applied this procedure for all the OBs, there should be a `shifts.txt` file in each of the `helpers` subdirectory, even if you kept the original shifts proposed by the pipeline. Now that each exposures have a correct relative astrometry within a given OB, we need to check the relative astrometry between the different OBs themselves. To do so, navigate back to the working directory and open all the combined images of your helper target in all OBs at once (don't forget to edit the `xxx`):
```bash
ds9 $(find | cat | sort | grep "helpers/combine_sci_reconstructed_xxx_img_cont.fits")
```

12) Repeat step (6) above to match all the frames.

13) In each frame, create a circular region centered on the peak of emission of your helper target. Open a new text file called `centroid_helper.txt` and note down the RA and Dec position of the center of this circle, in degrees, and in the right order (i.e., first frame first). If, for some reason, you cannot find your helper target in one of the OB, just write random coordinates: we will deal with this later on. The format is:
```bash
ra1 dec1
ra2 dec2
... ...
```

14) In the working directory, create a new directory called `sci-shift-master` and move the `centroid_helper.txt` file there. Also make a copy of `sci-master/reduce.sh` there and open this copy in your text editor. Change the first line to look like:
```bash
esorex kmos_combine -edge_nan -method='user' -filename='shifts.txt' combine.sof
```

15) Now copy the file `make_shifts.cpp` there and compile it:
```bash
cphy++ optimize make_shifts.cpp
```

16) This script will compile together all your previous shift measurements into a single master shift list that the pipeline can use. You call it this way:
```bash
./make_shifts "../sci-" helper=xxx
```
Make sure to replace `xxx` with the name of your helper target, and if needed change `"../sci-"` to match the starting pattern of the directories of your reduced OBs.

At this stage, you can also exclude one or several OBs from the reduction, in particular if you failed to detect your helper target or if the data quality is bad. To do so, simply specify their ID in the `exclude` argument. For example, if you want to exclude the OBs `sci-03` and `sci-10`, write instead:
```bash
./make_shifts "../sci-" helper=xxx exclude=[3,10]
```

This program will create the `shifts.txt` file and the SOF file for the pipeline. You can inspect these files if you wish, and then run the `reduce.sh` script. Note that this reduction may contain fewer objects than what you had in section (H). This can happen if a) you have put science targets inside the sky pointings, or b) if you have multiple target lists in your observing program. Indeed, when you provide manual position shifts, the pipeline will only reduce the exposures for which these shifts were computed (i.e., all the exposures containing your helper target). To reduce the missing targets, you will have to use another helper target that was observed simultaneously, and re-do all of this for these other exposures. If you have no other helper target, then there is nothing you can do but use the standard reduction without shifts for these objects.

17) Now you can repeat the step (H.5) to inspect the continuum images, and see if your shifts have improved the signal to noise ratio.

# Appendix A. Analyzing cubes: extracting spectra

To extract spectra from the cubes you can use QFitView (see section E), which only allows you to extract all the flux within a given pixel or circular aperture. The esorex pipeline and the recipe `kmos_extract_spec` is more flexible, as it allows you to use arbitrary masks. But you may want something more...

There are cases where you need to de-blend two close objects in the IFU and extract both their spectra. You may also want to take into account a uniform background level. To do this, you can use the `multispecfit.cpp` tool. Given a fixed set of 2D models, it will find the optimal linear combination of these models to describe the content of the KMOS data cube at each wavelength, independently. It will then give you the spectrum associated to each model.

1) First go into the directory where you have created your final, combined cubes in section (H). Create yourself a new directory called 'spectra'; go there and copy the `multispecfit.cpp` program and compile it:
```bash
cphy++ optimize multispecfit.cpp
```

2) Now you want to create the models that will describe your data. The simplest and most common case is to have a uniform background and a Gaussian profile. However, to demonstrate fully the usage of the tool, I will assume that your IFU contains two small galaxies close to one another. You can easily adapt this to your own situation.

3) Here is how you proceed to build the models with IDL (sorry Python folks):
```IDL
; Load the continuum image of your object 'xxx' just to get the right image dimensions.
; If you know 'nx' and 'ny', you can also specify them directly
img = mrdfits('../continuum/combine_sci_reconstructed_xxx_img_cont.fits')
nx = n_elements(img[*,0])
ny = n_elements(img[0,*])

; Create the 'x' and 'y' pixel coordinate arrays
pp = array_indices(img, findgen(nx, ny))
px = reform(pp[0,*], nx, ny)
py = reform(pp[1,*], nx, ny)

; We know (for example) that the first galaxy is centered on pixels (5,4), as seen
; in DS9 or QFitsView, while the second is centered on (7,6). Note that both DS9
; and QFitsView use 1-based pixel coordinates (the first pixel has coordinates 1,1),
; but here in IDL we need 0-based values, so we subtract 1 from the coordinates.
; Lastly, we assume that both galaxies are well described with a Gaussian of width 1.2
; pixels. It is up to you to figure out the best value, either by trying to match the
; continuum emission, or some line map you will have previously built, or even some
; ancillary imaging from another facility.

galaxy1 = exp(-((px - 4)^2 + (py - 3)^2)/(2.0*1.2^2))
galaxy2 = exp(-((px - 6)^2 + (py - 5)^2)/(2.0*1.2^2))

; NB: The models are re-normalized internally by the tool to make sure that they contain
; a flux of unity to start with. This means that the value of the resulting spectrum is
; indeed the total flux corresponding to the model.

; The last step is to store these models into a cube for the program
models = fltarr(3, nx, ny)
models[0,*,*] = 1       ; first a constant background
models[1,*,*] = galaxy1 ; the first galaxy
models[2,*,*] = galaxy2 ; the second galaxy

; And write that down
mwrfits, /create, models, 'combine_sci_reconstructed_xxx_models.fits'

```

4) Now that the models are created, we can call the program to do the fit:
```bash
./multispecfit ../combine_sci_reconstructed_xxx.fits combine_sci_reconstructed_xxx_models.fits suffix=[bg,gal1,gal2]
```

This will create three files:
```bash
combine_sci_reconstructed_xxx_bg_spec.fits
combine_sci_reconstructed_xxx_gal1_spec.fits
combine_sci_reconstructed_xxx_gal2_spec.fits
```
The first contains the background level, while the second and third contain the spectra of the two galaxies. The format of these FITS files is the same as what `kmos_extract_spec` would create: the first extension is empty, the second contains the spectrum and the third contains the uncertainty. The wavelength corresponding to each pixel is given by the WCS system.
