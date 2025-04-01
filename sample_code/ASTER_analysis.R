library(data.table)
library(meshed)

# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

# data description at https://lpdaac.usgs.gov/products/ag100v003/
# downloaded from https://search.earthdata.nasa.gov/search/granules?p=C2763266348-LPCLOUD&pg[0][v]=f&pg[0][gsk]=-start_date&g=G2817988042-LPCLOUD&sb[0]=-110.24121%2C40.80957%2C-108.5625%2C42.68214&fi=ASTER&fpj=ASTER%2BGED&tl=1731637591!3!!&lat=42.3544921875&long=-114.46875&zoom=6
load("data/ASTER.RData")

# NOTE FOR LATER: Phi prior Unif(0, 200)?

