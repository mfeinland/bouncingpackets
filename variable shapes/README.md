
# Variable Shapes Paper Summary
### Materials relevant to *Variable shapes observed in electron bouncing microburst packets* article, published 2024 in GRL

**Author:** Max Feinland for Blum Research Group, LASP \
**Contact:** max.feinland@lasp.colorado.edu or maxfeinland@gmail.com \
**Date Created:** 9/6/2024\
**Last Modified:** 9/6/2024


## Directory
- `all_events_v2.csv` has all 473 identified bouncing microburst packet candidates, with timestamp, lat/lon, MLT, L-shell, shape classification (if applicable) and quality classification.
- `bouncingpacketfunctions.py` has the functions (O'Brien algorithm, authored algorithm) necessary to search the SAMPEX data for bouncing packets.
- `findbouncingpackets.py` is the driver script that lets you query dates to search.
- `make_plots` has 4 functions to generate the 4 figures in the manuscript.
- `model_preds_and_spacings.csv` has model predictions and exact spacings of each "good" microburst, used to generate Figure 3 in the manuscript .
