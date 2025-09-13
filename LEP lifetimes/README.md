
# LEP Lifetimes Paper Summary
### Materials relevant to "Lightning-induced precipitation as a proxy for inner belt MeV electron decay times", submitted 2025 in *JGR Space Physics*

**Author:** Max Feinland for Blum Research Group, LASP \
**Contact:** max.feinland@lasp.colorado.edu or maxfeinland@gmail.com \
**Last Modified:** 9/12/2025


## Directory
- `find_bursts.py` searches through all of the SAMPEX/HILT State 4 data and applies a modified [O'Brien (2003)](https://doi.org/10.1029/2002JA009784) algorithm to detect bursts. There are then some quality checks that it applies. The final result is output to a `.csv` file containing the timestamps for all microbursts and ephemeris data. This is detailed in the `Data_Files` subfolder.
- `find_bouncing_packets.py` then restricts that catalog to *L* < 2.5 and plots each burst event, prompting the user to categorize into a bouncing packet or other burst event. The events identified as bouncing packets are saved to `inner_belt_bursts.csv`.
- `make_figures.py` has the code to replicate all the figures from the manuscript.

`Data_Files` has the following files:
  - `inner_belt_bursts.csv` contains the timestamps for all detected microbursts at *L* < 3 and the interpolated values for latitude, longitude, altitude, L-shell, local B field magnitude (`B`), MLT, losscone angle in the same hemisphere as the spacecraft (`losscone1`), losscone angle in either hemisphere (`losscone2`), HILT boresight pitch angle (`pitch`), the SAA flag, the attitude flag, and the count rate at that timestamp. 
  - `bouncing_packets.csv` has all the timestamps for the microbursts identified as bouncing packets and all the same ephemeris data as above, plus the equatorial B-field (`Beq`). 
  - `omni.csv` has the hourly Kp index, hourly Dst index, and sunspot number, taken from the [NASA OMNI database.](https://omniweb.gsfc.nasa.gov/)
  - `LISOTD_LRMTS_V2.3.2015.nc.zip` has monthly lightning data gridded in 2.5x2.5 degree cells from the LIS/OTD satellites. I zipped it so GitHub would accept it, but needs to be unzipped for data analysis. Downloaded from the [NASA Earthdata database.](https://www.earthdata.nasa.gov/data/catalog/ghrc-daac-lolrmts-2.3.2015)
  - `Lgrid.dat` has L-shell contours made by Sergio Vidal-Luengo using the [aacgmv2](https://pypi.org/project/aacgmv2/) Python package. I used them for Figure 4, which has a geographic map.
  - `ne_110m_admin_0_countries_lakes.shp` and `.shx` are the shapefiles to plot the world on a global map. This was taken from [Natural Earth](https://www.naturalearthdata.com/downloads/110m-cultural-vectors/).
  - `gaussian_fits.csv` has data describing the Gaussian fit of each consecutive bouncing packet peak that had *R*^2^ > 0.9: the timestamp (`t`), width (`w`) in seconds, spacing (`s`) in seconds, and error associated each spacing (`s_err`).
  - `missing_times.csv` and `spin.csv` are catalogs of epochs during which the spacecraft was spinning, the attitude data was flagged, or data were missing. This is to plot these epochs for the time series in Figure 5. 
