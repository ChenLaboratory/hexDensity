# hexDensity 1.4.10 [2025-09-28]

## Bug fixes and minor improvements:
* Better bound on numbers of columns. Previously, hexbin may add an extra column
out of bound to deal with points right on the edges. This is acceptable for 
binning but not great for KDE due to the need for edge correction.
* new NEWS.md.