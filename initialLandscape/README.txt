# .asc rasters
aspect.asc       (degrees)
cellID.asc       raster cell ID
elev.asc         elevation (m a.s.l)
forestMask.asc   forest extent within the park
                 0 = outside forest
                 1 = inside forest
parkMask.asc     park extent
                 0 = outside park
                 1 = inside park
slope.asc        (degrees)
swhc.asc         soil water holding capacity (cm)



# trees75.csv
cellID:      raster cell ID
sp:          species latin name
n:           number of trees
dbh:         diameter at breast height (cm)



# envVariables.csv
cellID:      raster cell ID
park:        park extent
             0 = outside park
             1 = inside park
forest:      forest extent within the park
             0 = outside forest
             1 = inside forest
elev:        elevation (m a.s.l.)
slope:       (degrees)
aspect:      (degrees)
swhc:        soil water holding capacity (cm)
GRECO:       GRECO region ID (defined by French NFI)
Cd_crbn:     carbonate mineral content of the bedrock
             0 = no or few carbonates
             1 = moderate presence of carbonates
             2 = rich in carbonates
             3 = Nearly pure carbonates
Cd_hydr:     permeability potential of the bedrock according to its mineralogical nature
             0 = Highly permeable
             1 = Permeable
             2 = Low permeability
             3 = Impermeable or nearly impermeable
SIQpet:      salem SI (site index) prediction for Q. petraea
SIFsyl:      salem SI (site index) prediction for F. sylvatica
SIAalb:      salem SI (site index) prediction for A. alba
SIPabi:      salem SI (site index) prediction for P. abies
