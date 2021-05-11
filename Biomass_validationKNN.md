---
title: "Biomass_validationKNN"
author: "Ceres Barros and Alex Chubaty"
date: "11 May 2021"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---



# Overview

Provide an overview of what the module does / how to use the module.

Module documentation should be written so that others can use your module.
This is a template for module documentation, and should be changed to reflect your module.

## R Markdown

R Markdown syntax allows R code, outputs, and figures to be rendered in the documentation.

For help writing in R Markdown, see https://rmarkdown.rstudio.com/.

# Usage


```r
library(SpaDES.core)

setPaths(modulePath = file.path(".."))
getPaths() # shows where the 4 relevant paths are

times <- list(start = 0, end = 10)

parameters <- list(
  #.progress = list(type = "text", interval = 1), # for a progress bar
  ## If there are further modules, each can have its own set of parameters:
  #module1 = list(param1 = value1, param2 = value2),
  #module2 = list(param1 = value1, param2 = value2)
)
modules <- list("Biomass_validationKNN")
objects <- list()
inputs <- list()
outputs <- list()

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects)

mySimOut <- spades(mySim)
```

# Parameters

Provide a summary of user-visible parameters.


```
## defineParameter: 'coverThresh' is not of specified type 'integer'.
```

```
## defineParameter: 'LCChangeYr' is not of specified type 'integer'.
```

```
## defineParameter: 'validationYears' is not of specified type 'integer'.
```

```
## defineParameter: '.useCache' is not of specified type 'logical'.
```

```
## Loading required package: dplyr
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```



|paramName              |paramClass |default      |min  |max  |paramDesc                                                                                                                                                                                                                                                                      |
|:----------------------|:----------|:------------|:----|:----|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|coverThresh            |integer    |10           |NA   |NA   |The minimum % cover a species needs to have (per pixel) in the study area to be considered present. Should be the same as the one used to obtain the species cover layers for simulation set up.                                                                               |
|deciduousCoverDiscount |numeric    |0.8418911    |NA   |NA   |This was estimated with data from NWT on March 18, 2020 and may or may not be universal. Should be the same as the one used when preparing 'cohortData' in the simulation set up.                                                                                              |
|LCChangeYr             |integer    |2001, 20.... |1985 |2015 |An integer or vector of integers of the validation period years, defining which years of land-cover changes (i.e. disturbances) should be excluded. Only used if rstLCChangeYr is not NULL. See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information. |
|minCoverThreshold      |numeric    |5            |0    |100  |Cover that is equal to or below this number will be omitted from the dataset Should be the same as the one used when preparing 'cohortData' in the simulation set up.                                                                                                          |
|obsDeltaAgeB           |logical    |TRUE         |NA   |NA   |When TRUE, the observed changes in biomass and age (deltaB, deltaAge) between the two validation years will be plotted as maps and scatterplots                                                                                                                                |
|pixelGroupBiomassClass |numeric    |100          |NA   |NA   |When assigning pixelGroup membership, this defines the resolution of biomass that will be considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same Should be the same as the one used when preparing 'cohortData' in the simulation set up.   |
|sppEquivCol            |character  |Boreal       |NA   |NA   |The column in sim$specieEquivalency data.table to use as a naming convention                                                                                                                                                                                                   |
|validationReps         |integer    |1, 2, 3,.... |NA   |NA   |The simulation repetitions for the validation. Defaults to 1:10. Set to NA if not using repetitions (i.e. only one run)                                                                                                                                                        |
|validationYears        |integer    |2001, 2011   |NA   |NA   |The simulation years for the validation. Defaults to 2001 and 2011. Must select two years                                                                                                                                                                                      |
|.plotInitialTime       |numeric    |0            |NA   |NA   |Vector of length = 1, describing the simulation time at which the first plot event should occur. Set to NA to turn plotting off.                                                                                                                                               |
|.plotInterval          |numeric    |NA           |NA   |NA   |This describes the simulation time interval between plot events                                                                                                                                                                                                                |
|.saveInitialTime       |numeric    |NA           |NA   |NA   |This describes the simulation time at which the first save event should occur                                                                                                                                                                                                  |
|.savePlots             |logical    |TRUE         |NA   |NA   |Whether plots should be saved in file.path(outputPath(sim), 'Figs')                                                                                                                                                                                                            |
|.saveInterval          |numeric    |NA           |NA   |NA   |This describes the simulation time interval between save events                                                                                                                                                                                                                |
|.studyAreaName         |character  |NA           |NA   |NA   |Human-readable name for the study area used. If NA, a hash of studyArea will be used.                                                                                                                                                                                          |
|.useCache              |logical    |init         |NA   |NA   |Controls cache; caches the init event by default                                                                                                                                                                                                                               |

# Events

Describe what happens for each event type.

## Plotting

Write what is plotted.

## Saving

Write what is saved.

# Data dependencies

## Input data

How to obtain input data, and a description of the data required by the module.
If `sourceURL` is specified, `downloadData("Biomass_validationKNN", "..")` may be sufficient.


```
## defineParameter: 'coverThresh' is not of specified type 'integer'.
```

```
## defineParameter: 'LCChangeYr' is not of specified type 'integer'.
```

```
## defineParameter: 'validationYears' is not of specified type 'integer'.
```

```
## defineParameter: '.useCache' is not of specified type 'logical'.
```



|objectName         |objectClass              |desc                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |sourceURL                                                                                                                                                                                                      |
|:------------------|:------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|allCohortData      |data.table               |All 'cohortData' tables saved during the simulation, particularly for the validation years. If not supplied, the module will attempt to retrieve them using the 'simulationOutputs' table                                                                                                                                                                                                                                                                                        |NA                                                                                                                                                                                                             |
|biomassMap         |RasterLayer              |total biomass raster layer in study area, filtered for pixels covered by cohortData. Only used to calculate total no. of pixels being simulated If not supplied, will default to to the Canadian Forestry Service, National Forest Inventory, kNN-derived total aboveground biomass map from 2001. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata                                                                                  |NA                                                                                                                                                                                                             |
|firePerimeters     |sf                       |A map of fire perimeters in the study area that can be used to exclude pixels that have been burnt during the validation period. Defaults to the Canadian Wildland Fire Information System 1986-2018 National Burned Area Composite, subset to fires between 2001 and 2011 (inclusively).                                                                                                                                                                                        |http://cwfis.cfs.nrcan.gc.ca/downloads/nbac/nbac_1986_to_2019_20200921.zip                                                                                                                                     |
|pixelGroupMapStk   |RasterStack              |A stack of pixelGroupMaps saved during the simulation, particularly for the validation years. If not supplied, the module will attempt to make it using the 'simulationOutputs' table                                                                                                                                                                                                                                                                                            |NA                                                                                                                                                                                                             |
|rawBiomassMapStart |RasterLayer              |observed total biomass raster layer in study area at the first year of the validation period. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived total aboveground biomass map from 2001 (in tonnes/ha). If necessary, biomass values are rescaled to match changes in resolution. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata.                                                                  |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif |
|rawBiomassMapEnd   |RasterLayer              |observed total biomass raster layer in study area at the last year of the validation period. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived total aboveground biomass map from 2011 (in tonnes/ha). If necessary, biomass values are rescaled to match changes in resolution. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990                                                                                 |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif |
|rasterToMatch      |RasterLayer              |A raster of the studyArea in the same resolution and projection as rawBiomassMapStart. This is the scale used for all outputs for use in the simulation.                                                                                                                                                                                                                                                                                                                         |NA                                                                                                                                                                                                             |
|rstLCChange        |RasterLayer              |A mask-type map of land cover changes in the study area that can be used to exclude pixels that have been disturbed during the validation period. Defaults to Canada's forest change map between 1985-2011 (CFS), filtered for years 2001-2011 (inclusively) and all disturbances collapsed (map only has values of 1 and NA). See parameter LCChangeYr to change the period of disturbances, and https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information. |https://opendata.nfis.org/downloads/forest_change/C2C_change_type_1985_2011.zip                                                                                                                                |
|rstLCChangeYr      |RasterLayer              |An OPTIONAL map of land cover change years in the study area used to exclude pixels that have been disturbed during the validation period. Defaults to Canada's forest change national map between 1985-2011 (CFS). By default disturbances are subset to to years 2001-2011 (inclusively; see parameter LCChangeYr). See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information.                                                                         |https://opendata.nfis.org/downloads/forest_change/C2C_change_year_1985_2011.zip                                                                                                                                |
|simulationOutputs  |data.table               |An OPTIONAL table listing simulation outputs (as passed to 'spades()', or 'experiment that will be used to make 'allCohortData', 'pixelGroupMapStk', if these are not provided.                                                                                                                                                                                                                                                                                                  |NA                                                                                                                                                                                                             |
|speciesLayersStart |RasterStack              |observed cover percentage raster layers by species in Canada species map, at the first year of the validation period. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived species cover maps from 2001, using a cover threshold of 10% - see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata                                                                                                              |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/                                                                     |
|speciesLayersEnd   |RasterStack              |observed cover percentage raster layers by species in Canada species map used for validation at the last year of the validation period. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived species cover maps from 2011 - see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata                                                                                                                            |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/                                                                     |
|sppColorVect       |character                |named character vector of hex colour codes corresponding to each species                                                                                                                                                                                                                                                                                                                                                                                                         |NA                                                                                                                                                                                                             |
|sppEquiv           |data.table               |table of species equivalencies. See LandR::sppEquivalencies_CA.                                                                                                                                                                                                                                                                                                                                                                                                                  |NA                                                                                                                                                                                                             |
|standAgeMapStart   |RasterLayer              |observed stand age map in study area, at the first year of the validation period Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived biomass map from 2001 - see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata                                                                                                                                                                                          |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif                    |
|standAgeMapEnd     |RasterLayer              |observed stand age raster layer in study area, at the last year of the validation period. Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived stand age map from 2011. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990                                                                                                                                                                                             |http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/NFI_MODIS250m_2011_kNN_Structure_Stand_Age_v1.tif                    |
|studyArea          |SpatialPolygonsDataFrame |Polygon to use as the study area. Defaults to an area in Southwestern Alberta, Canada.                                                                                                                                                                                                                                                                                                                                                                                           |NA                                                                                                                                                                                                             |

## Output data

Description of the module outputs.


```
## defineParameter: 'coverThresh' is not of specified type 'integer'.
```

```
## defineParameter: 'LCChangeYr' is not of specified type 'integer'.
```

```
## defineParameter: 'validationYears' is not of specified type 'integer'.
```

```
## defineParameter: '.useCache' is not of specified type 'logical'.
```



|objectName         |objectClass |desc                                                                                                                                                                                              |
|:------------------|:-----------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|rstDisturbedPix    |RasterLayer |Raster of pixel IDs (as a mask) that have been disturbed by fire or suffered land-cover changes during the validation period. These pixels are excluded form the validation.                      |
|rawBiomassMapStart |RasterLayer |observed total biomass raster layer in study area at the first year of the validation period. Filtered to exclude pixels that were disturbed during the validation period                         |
|rawBiomassMapEnd   |RasterLayer |observed total biomass raster layer in study area at the last year of the validation period. Filtered to exclude pixels that were disturbed during the validation period                          |
|speciesLayersStart |RasterStack |observed cover percentage raster layers by species in Canada species map, at the first year of the validation period. Filtered to exclude pixels that were disturbed during the validation period |
|speciesLayersEnd   |RasterStack |observed cover percentage raster layers by species in Canada species map, at the last year of the validation period. Filtered to exclude pixels that were disturbed during the validation period  |
|standAgeMapStart   |RasterLayer |observed stand age map in study area, at the first year of the validation period Filtered to exclude pixels that were disturbed during the validation period                                      |
|standAgeMapEnd     |RasterLayer |observed stand age map in study area, at the last year of the validation period Filtered to exclude pixels that were disturbed during the validation period                                       |
|standCohortData    |data.table  |                                                                                                                                                                                                  |

# Links to other modules

Describe any anticipated linkages to other modules.
