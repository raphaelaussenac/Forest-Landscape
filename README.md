# Goal
Assign trees (specifying their species, diameter and height) to all 25m*25m forest cells over different landscapes from LIDAR and inventory data.



# Outputs

## rasters
* aspect.asc: aspect (degrees)
* cellID25.asc: 25m*25m cell ID
* cellID100.asc: 100m*100m cell ID
* elev.asc: elevation (m a.s.l)
* parkMask.asc: Bauges Geopark extent
  + 0 = outside park
  + 1 = inside park
* pH.ascii: soil pH
* slope.asc: slope (degrees)
* swhc.asc: soil water holding capacity (cm)


## trees.csv
* cellID25: 25m*25m cell ID
* sp: species latin name
* n: number of trees
* dbh: diameter at breast height (cm)
* h: height (m)


## cell25.csv
* cellID25: 25m*25m cell ID
* cellID100: 100m*100m cell ID
* park: Bauges Geopark extent
  + 0 = outside park
  + 1 = inside park
* elev: elevation (m a.s.l.)
* slope: slope (degrees)
* aspect: aspect (degrees)
* swhc: soil water holding capacity (cm)
* pH: soil pH
* GRECO: GRECO region ID (defined by French NFI)
* SIQpet: salem SI (site index) prediction for Q. petraea
* SIFsyl: salem SI (site index) prediction for F. sylvatica
* SIAalb: salem SI (site index) prediction for A. alba
* SIPabi: salem SI (site index) prediction for P. abies


## managTableCell100.csv
* cellID100: 100m*100m cell ID
* owner: ownership (public or private)
* access: accessibility
  + 0 = inaccessible area
  + 1 = accessible area
* protect: protection status
  + 0 = non-protected area
  + 1 = protect area
* forestCellsPerHa: nb of 25m*25m cells with forest within 100m *100m cells
* compoType: composition type
* gini: gini index (calculated on individual tree basal area)
* rdi: rdi density index
* BA: total basal area (m2)
* BA_ha: basal area per (m2/ha). May be different from 'BA' depending on 'forestCellsPerHa'
* Dg: mean quadratic diameter at breast height (cm)
* meanH: mean height (m)
* structure: even-aged or uneven-aged structure
* density: density level (high, medium, low)
* manag: stand management (management - density)
