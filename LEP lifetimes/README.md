# LEP Lifetimes Paper Summary
### Materials relevant to *"Lightning-induced relativistic electron precipitation as a proxy for inner belt MeV electron lifetimes"*, to be submitted 2025 in *JGR Space Physics*

**Author:** Max Feinland for Blum Research Group, LASP \
**Contact:** max.feinland@lasp.colorado.edu or maxfeinland@gmail.com \
**Last Modified:** 5/19/2025


## Directory
- `find_microbursts.py` searches through all of the SAMPEX/HILT State 4 data and applies a modified [O'Brien (2003)](https://doi.org/10.1029/2002JA009784) algorithm to detect microbursts. There are then some quality checks that it applies. The final result is output to a .csv file containing the timestamps for all microbursts and ephemeris data. This is detailed in the `Data_Files` subfolder.
- `Data_Files` has:
      -    `inner_belt_microbursts.csv`, which contains the timestamps for all microbursts and the interpolated values for latitude, longitude, altitude, L-shell, local B field magnitude (`B`), MLT, losscone angle in the same hemisphere as the spacecraft (`losscone1`), losscone angle in either hemisphere (`losscone2`), HILT boresight pitch angle (`pitch`), the SAA flag, the attitude flag, and the count rate at that timestamp. 
