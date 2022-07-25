# Goal
Assign trees (specifying their species, diameter and height) to all 25m*25m forest cells over different landscapes from LIDAR and inventory data.



# Outputs

## rasters
* aspect.asc: aspect (degrees)
* cellID25.asc: 25m*25m cell ID
* cellID100.asc: 100m*100m cell ID
* elev.asc: elevation (m a.s.l)
* parkMask.asc: case study area extent
  + 0 = outside case study area
  + 1 = inside case study area
* forestMask: forest cells to be simulated
  + 0 = not simulated (no trees at the initial state or outside the simulated area)
  + 1 = simulated (presence of trees at the initial state and inside the simulated area)
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
* park: case study area extent
  + 0 = outside case study area
  + 1 = inside case study area
* forest: forest cells to be simulated
  + 0 = not simulated (no trees at the initial state or outside the simulated area)
  + 1 = simulated (presence of trees at the initial state and inside the simulated area)
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
* dgModul: values (-5/0/+5) in cm to be added to the prescribed Dg (used for the 'working for complexity' alternative management)

Alternative landscapes are created to simulate the alternative managements. Letters following the file name correspond to the alternative managements:
* B: baseline
* BC: baseline + working for complexity
* E: extensification
* EC: extensification + working for complexity
* I: intensification
* IC: intensification + working for complexity

Differences between the baseline landscapes and the alternative landscapes:

Bauges:
* intensification: -8% of area to unmanaged<sup>*</sup>.
* extensification: +8% of area to unmanaged<sup>*</sup>.
* working for complexity:
  + increase silvicultural system evenness in the landscape by increasing the Gini threshold defining the limit between even- and uneven-aged stands untill reaching a balanced landscape (in term of even- vs uneven-aged stands).
  + modulate prescribed Dg/Dharv values by adding -5, 0 or +5 cm to the prescribed Dg/Dharv values evenly across composition types (use dgModul column for that).

<sup>*</sup> In the Bauges Landscape, one part of the forest cannot be managed (it is either not accessible and or on too steep slopes), the other part recieves a value of skidding distance (SkD). In the baseline landscape, stands with an associated SkD <= 1500m can be managed. In the intensification landscape, all stands with a SkD can be managed, which amounts to a reduction of the surface of unmanaged forests of about 8%. In the extensification landscape, we decreased the SkD threshold below which forests can be managed so as to get an increase of the surface of unmanaged forests of about 8% (SkD = 1000).

Milicz:
* intensification: 0% of area to unmanaged.
* extensification: 20% of area to unmanaged.
* working for complexity:
  + manage quercus and fagus stands as uneven-aged stands.
  + modulate prescribed Dg/Dharv values by adding -5, 0 or +5 cm to the prescribed Dg/Dharv values evenly across composition types (use dgModul column for that).

Sneznik:
* intensification: 0% of area to unmanaged.
* extensification: 20% of area to unmanaged.
* working for complexity:
  + increase silvicultural system evenness in the landscape by decreasing the Gini threshold defining the limit between even- and uneven-aged stands untill reaching a balanced landscape (in term of even- vs uneven-aged stands).
  + modulate prescribed Dg/Dharv values by adding -5, 0 or +5 cm to the prescribed Dg/Dharv values evenly across composition types (use dgModul column for that).