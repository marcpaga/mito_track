# Mito tracking

This repository contains the code used for the analysis of the mitochondrial flow in CITATION NEEDED.

There are four scripts:
- `mito_track.py`: corner detection and flow estimation using the OpenCV library. It outputs several .csv files with the coordinates of all the tracked points.
- `coord_process.R`: estimates the membrane location, direction and distance to the membrane of the tracked points.
- `plot.R`: code necessary to produce the same plots as in Supplemental Figure N.
- `Mito_to_membrane_distance.m`:calculates the distance between mitochondria and membrane over 360 degrees.

# Dependencies

- `Python` = `3.7.2`
    - `open-cv` = `4.1.2.30`
    - `numpy` = `1.17.4`
    - `tiffcapture` = `0.1.6`
- `R` = `3.6.1`
    - `ijtiff` = `2.0.4`
    - `ggplot2` = `3.2.1`
    - `RColorBrewer` = `1.1`
    - `tiff` = `0.1`
- `Matlab` = `Mathworks version R2014a`

# License

[Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)