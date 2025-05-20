
# LEP Lifetimes Paper Summary
### Materials relevant to "Lightning-induced precipitation as a proxy for inner belt MeV electron lifetimes", to be submitted 2025 in *JGR Space Physics*

**Author:** Max Feinland for Blum Research Group, LASP \
**Contact:** max.feinland@lasp.colorado.edu or maxfeinland@gmail.com \
**Last Modified:** 5/19/2025


## Directory
- `find_microbursts.py` searches through all of the SAMPEX/HILT State 4 data and applies a modified [O'Brien (2003)](https://doi.org/10.1029/2002JA009784) algorithm to detect microbursts. There are then some quality checks that it applies. The final result is output to a .csv file containing the timestamps for all microbursts and ephemeris data. This is detailed in the `Data_Files` subfolder.
- `find_bouncing_packets.py` does **x, y, and z.**
- `Data_Files` has the following files:
  - `inner_belt_microbursts.csv` contains the timestamps for all microbursts and the interpolated values for latitude, longitude, altitude, L-shell, local B field magnitude (`B`), MLT, losscone angle in the same hemisphere as the spacecraft (`losscone1`), losscone angle in either hemisphere (`losscone2`), HILT boresight pitch angle (`pitch`), the SAA flag, the attitude flag, and the count rate at that timestamp. 
  - `bouncing_packets.csv` has all the timestamps for the microbursts identified as bouncing packets and all the same ephemeris data as above, plus the equatorial B-field (`Beq`). 
  - `omni.csv` has the hourly Kp index, hourly Dst index, and sunspot number **TAKEN FROM??**
  - `LISOTD_LRMTS_V2.3.2015.nc.zip` has monthly lightning data gridded in 2.5x2.5 degree cells from the LIS/OTD satellites. I zipped it so GitHub would like it, but needs to be unzipped for data analysis. **TAKEN FROM??**
  - `Lgrid.dat` has L-shell contours made by Sergio Vidal-Luengo using the [aacgmv2](https://pypi.org/project/aacgmv2/) Python package. I used them for Figure 4, which has a geographic map.
  - Finally, `ne_110m_admin_0_countries_lakes.shp` has the shapefile to plot the world on a global map. This was taken from Natural Earth and is available for download [here.](https://www.naturalearthdata.com/downloads/110m-cultural-vectors/)
