# LEP Paper Summary
### Materials relevant to *Lightning-induced relativistic electron precipitation from the inner radiation belt* article, published 2024 in *Nature Communications*

**Author:** Max Feinland for Blum Research Group, LASP \
**Contact:** max.feinland@lasp.colorado.edu or maxfeinland@gmail.com \
**Last Modified:** 9/3/2024


## Directory
- **Drafts & Reviewer Comments** has just that; there are 4 paper versions for 3 rounds of review, so paper draft version 2 was prompted by round 1 of review.
- **Generate Plots** has materials needed to replicate plots from the paper; it's a MATLAB script that reads in the data from the excel file.
- **LaTeX** has the materials needed to recompile the paper in Overleaf or another LaTeX processing software/service/etc. It has a .tex document (the main manuscript file), all the figures referenced in the manuscript file, style file so it compiles how Nature wants and a trackchanges file so the `\add`, `\remove`, `\change`, and `\note` commands (which are part of the `trackchanges` package) don't throw errors.
- `Low_L_Driver_Script.m` has the code necessary to search through SAMPEX data and identify bouncing packet microbursts. If you were to run this from August 7th 1996 to August 7th 2006, you should find 68 events. Of those, 45 events are likely bouncing microburst packets.
    - This code downloads HILT count rate data and SAMPEX ephemeris data from online. It auto-deletes these files if it needed to download them, but you can change this if you want in the code. ==add option in code for this==
- `lshelldata.csv` has these events and their properties (timestamp of  microburst beginning, average period, latitude, longitude, L-shell, magnetic local time, maximum count rate, SAA flag (as determined by SAMPEX ephemeris data) )
