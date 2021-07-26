# Goal
Assign trees (specifying their species and size) to all 25*25m forest cells over different landscapes from LIDAR and inventory data.



# outputs

## rasters
* aspect.asc       (degrees)
* cellID25.asc     25m*25m cell ID
* cellID100.asc    100m*100m cell ID
* elev.asc         elevation (m a.s.l)
* forestMask.asc   forest extent within the park
                   + 0 = outside forest
                   + 1 = inside forest
* parkMask.asc     park extent
                   + 0 = outside park
                   + 1 = inside park
* pH.ascii         soil pH
* slope.asc        (degrees)
* swhc.asc         soil water holding capacity (cm)


## trees75.csv
* cellID25:    25m*25m cell ID
* cellID100:   100m*100m cell ID
* sp:          species latin name
* n:           number of trees
* dbh:         diameter at breast height (cm)
* h:           height (m)


## envVariables.csv
* cellID25:           25m*25m cell ID
* cellID100:          100m*100m cell ID
* park:               park extent
                      + 0 = outside park
                      + 1 = inside park
* forest:             forest extent within the park
                      + 0 = outside forest
                      + 1 = inside forest
* elev:               elevation (m a.s.l.)
* slope:              (degrees)
* aspect:             (degrees)
* swhc:               soil water holding capacity (cm)
* pH:                 soil pH
* GRECO:              GRECO region ID (defined by French NFI)
* SIQpet:             salem SI (site index) prediction for Q. petraea
* SIFsyl:             salem SI (site index) prediction for F. sylvatica
* SIAalb:             salem SI (site index) prediction for A. alba
* SIPabi:             salem SI (site index) prediction for P. abies
* forestCellsPerHa:   nb of 25*25m cells with forest within 100*100m cells


## managTable.csv
* cellID100:          100m*100m cell ID
* protect:            protection status
                      + 0 = non-protected area
                      + 1 = protect area
* gini:               gini index (calculated on individual tree basal area)
* BA:                 total basal area (m2)
* Dg:                 mean quadratic diameter at breast height (cm)
* meanH:              mean height (m)
* structure:          even-aged or uneven-aged structure
* compoType:          composition type
* owner:              ownership (public or private)
* access:             accessibility
                      + 0 = inaccessible area
                      + 1 = accessible area
* forestCellsPerHa:   nb of 25*25m cells with forest within 100*100m cells
* rdi:                rdi density index
* density:            density level (high, medium, low)
* stand:              stand full management name
