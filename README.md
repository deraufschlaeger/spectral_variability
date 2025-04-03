# spectral_variability

In this repository can the code be found used in the bachelor's thesis of Birger Aufschläger.

My workflow for every object was as follows:
1. reduce spectrum with IRAF code
2. read out_obj.fits with the Blaze.ipynb
3. Check *every* order and modify Blaze function accordingly
4. run merge_orders function from wichtige_funktionen.py file with copy and pasted modified Blaze_estimate_bins function.
    This produces a .txt file with wavelength in the first column and corresponding weighted flux in the second

Diverting workflow for nova and other objects
Nova:
5. Use analysis_nova.ipynb notebook for analysis, plotting and making the latex table
other objects:
5. Use analysis_object.ipynb notebook for analysis (here example from the object TYC4454-1229-1), plotting an producing the latex table
