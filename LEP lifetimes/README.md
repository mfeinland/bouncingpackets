# LEP Lifetimes Paper Summary
### Materials relevant to *Lightning-induced relativistic electron precipitation as a proxy for inner belt MeV electron lifetimes* article, to be submitted 2025 in *JGR Space Physics*

**Author:** Max Feinland for Blum Research Group, LASP \
**Contact:** max.feinland@lasp.colorado.edu or maxfeinland@gmail.com \
**Last Modified:** 5/19/2025


## Directory
- `Generate_Plots` replicates plots from the paper; it's a MATLAB script that reads in the data from the excel file `Source_Data.xlsx.`
- `Low_L_Driver_Script.m` has the code necessary to search through SAMPEX data and identify bouncing packet microbursts. If you were to run this from August 7th 1996 to August 7th 2006, you should find 68 events. Of those, 45 events are likely bouncing microburst packets.
    - This code downloads HILT count rate data and SAMPEX ephemeris data from online, runs the [O'Brien (2003)](https://doi.org/10.1029/2002JA009784) algorithm to detect microbursts, and then runs an authored algorithm to find bouncing packets meeting certain criteria. There's a few options in there to make searching easier. Please report any issues with the code and I will try to fix them.
- `lshelldata.csv` has these events and their properties (timestamp of  microburst beginning, average period, latitude, longitude, L-shell, magnetic local time, maximum count rate, SAA flag (as determined by SAMPEX ephemeris data) )
