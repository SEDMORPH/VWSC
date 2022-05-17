# VWSC
Supercolour code for UDS and COSMOS

Please cite the following papers if you use this code:
1) Wild et al. 2014 http://adsabs.harvard.edu/abs/2014MNRAS.440.1880W
2) If you use the UDS code : Wild et al. 2016 http://adsabs.harvard.edu/abs/2016MNRAS.463..832W
3) If you use the COSMOS code : Wilkinson et al. 2021 https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.4533W

Start with the WRAPPER file and pay attention to the keywords to switch OFF different aspects of the code. You will need to edit the directories within this wrapper.

The code to create your own supercolour basis (eigenvectors) is included in this package, however the "stochastic burst" BC03 models which are used to do this are not available here. You can of course use any reasonable (or not so reasonable) set of stellar population models that vaguely cover galaxy colour parameter space. Some editing of the code will be required, in particular vwsc_read_stoch.pro . Please contact vw8@st-andrews.ac.uk if you would like help with this. 

If you are happy to use the eigenbases created by us, for the UDS and COSMOS, these are included. BE VERY CAREFUL TO ENSURE YOU EXACTLY MATCH THE FILTERS! 





