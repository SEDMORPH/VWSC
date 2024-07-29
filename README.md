# VWSC
Supercolour code for UDS and COSMOS

Please cite the following papers if you use this code:
1) Wild et al. 2014 http://adsabs.harvard.edu/abs/2014MNRAS.440.1880W
2) If you use the UDS code : Wild et al. 2016 http://adsabs.harvard.edu/abs/2016MNRAS.463..832W
3) If you use the COSMOS code : Wilkinson et al. 2021 https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.4533W

Start with the WRAPPER file and pay attention to the keywords to switch OFF different aspects of the code. You will need to edit the directories within this wrapper.

The code to create your own supercolour basis (eigenvectors) is included in this package, however the "stochastic burst" BC03 models which are used to do this are not available here. You can of course use any reasonable (or not so reasonable) set of stellar population models that vaguely cover galaxy colour parameter space. Some editing of the code will be required, in particular vwsc_read_stoch.pro . Please contact vw8@st-andrews.ac.uk if you would like help with this. 

If you are happy to use the eigenbases created for the papers above (for the UDS and COSMOS), these are available at [http://www-star.st-and.ac.uk/~vw8/downloads/](http://star-www.st-andrews.ac.uk/~web2vw8/downloads/index.html). BE VERY CAREFUL TO ENSURE YOU EXACTLY MATCH THE FILTERS IN YOUR DATASET TO THOSE USED HERE! 

Code dependencies:
* You will also need to install the VWPCA package that can be found in the same SEDMORPH GitHub repository.
* You will need the IDL coyote library, which can be found here: http://www.idlcoyote.com/documents/programs.php
* You will need the TEXtoIDL code http://physics.mnstate.edu/craig/textoidl/
* You will need to set you IDL_PATH environment variable, so that IDL knows where to look for all the code you install.






