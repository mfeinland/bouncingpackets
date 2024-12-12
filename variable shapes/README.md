# Variable Shapes Paper Summary
### Materials relevant to *Variable shapes observed in electron bouncing microburst packets* article, published 2024 in GRL

**Author:** Max Feinland for Blum Research Group, LASP \
**Contact:** max.feinland@colorado.edu or maxfeinland@gmail.com \
**Date Created:** 9/6/2024\
**Last Modified:** 12/12/2024

## Directory
- `bouncingpacketfunctions.py` has the functions (O'Brien algorithm, authored algorithm) necessary to search the SAMPEX data for bouncing packets.
- `findbouncingpackets.py` is the driver script that lets you query dates to search.
- `make_plots.py` has functions to generate the 4 figures in the manuscript.
- In the `Data_Files` folder:
    - `good_events.csv` has the 121 good bouncing microburst packets, with timestamp, mean period, predicted period, L-shell, MLT, lat, lon, shape category, quality score (1-3), portion of HILT field-of-view in the losscone for particles in the same hemisphere, spacings ($\Delta t$) for that event, and magnetic field model estimates (Tsyganenko-Sitnov 2005, Olson & Pfitzer quiet, Schulz & Lanzerotti dipole, Tsyganenko 1989, Ostapenko & Maltsev) for that event.
    - `all_events_v3.csv` has all 473 identified bouncing microburst packet candidates, with timestamp, lat/lon, MLT, L-shell, shape classification (if applicable) and quality classification.
    - `microburst_catalog_00.txt.zip` has Mike Shumko's microburst catalog made from the SAMPEX data. I had to zip it because GitHub didn't like the file size. You can also download that file from Mike's GitHub directly [here.](https://github.com/mshumko/sampex_microburst_widths/blob/main/sampex_microburst_widths/data/microburst_catalog_00.csv) This is used in Figure 3 (MLT/L distribution).
    - `Lgrid.dat` has L-shell contours made by Sergio Vidal-Luengo using the [aacgmv2](https://pypi.org/project/aacgmv2/) Python package. I used them for Figure 4, which has a geographic map.
    - Finally, `ne_110m_admin_0_countries_lakes.shp` has the shapefile to plot the world on a global map. This was taken from Natural Earth and is available for download [here.](https://www.naturalearthdata.com/downloads/110m-cultural-vectors/)
