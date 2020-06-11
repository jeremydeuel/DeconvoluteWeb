# SpectralDeconvolutionWeb
Javascript Application that allows the user to deconvolute UV-Vis spectra using the NNLS algorithm. Accepts spectra exported from an Agilent Cary 60 or a Thermo Nanodrop Spectrophotometer
# Install
Install all the files in this repository on a file accessible from the web. It is not possible to run this application locally, since most browsers reasonably deny the loading of local files.
# Configuration
The reference spectra to be used as well as the the colors associated with the substances are defined in the file refspectra.json
Each reference spectrum has a name and is defined by a dictionary with the following entries
* url: URL to the reference spectrum in csv format (see below)
* selected: Boolean, defines if this reference spectrum is loaded by default
* color: Default color for the graphical representation
# Reference Spectra
Reference Spectra should be recorded between 250 and 800nm with 1nm Step and 0.1s integration time in PBS pH 7.45 at a concentration of 6-7µM (maximum absorption 1.0) and at 37˚C. 
The Spectra are then exported (Agilent Cary 60) as CSV and the file edited in a text editor. The second line of the file has to indicate the concentration in µM by stating Conc: XXX µM (XXX being the numeric concentration, eg. 6.500)
# Dependencies
Dependencies are installed in the external folder
* Bootstrap v4.4.1 https://www.getbootstrap.com
* Numeric.js https://github.com/sloisel/numeric
* jQuery v3.4.1 http://www.jquery.com
* Chartjs v.2.6.0 https://github.com/chartjs/Chart.js
# Licence
GPL 3