# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_validationKNN",
  description = "Validation module for LandR Biomass predictions of forest succession. Based on Canadian Forest Service KNN maps", #"insert module description here",
  keywords = c("validation", "ecological simulation model",
               "forest dynamics", "forest succession", "data", "prediction"),
  authors = c(person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(Biomass_validationKNN = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_validationKNN.Rmd"),
  reqdPkgs = list("crayon", "raster", "achubaty/amc",
                  "sf", "XML", "RCurl", "ggplot2", "ggpubr",
                  "PredictiveEcology/LandR@modelBiomass (>=0.0.12.9003)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@development (>= 1.2.6.9005)",
                  "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development"),
  parameters = rbind(
    defineParameter("coverThresh", "integer", "10", NA, NA,
                    desc = paste("The minimum % cover a species needs to have (per pixel) in the study",
                                 "area to be considered present. Should be the same as the one used to obtain",
                                 "the species cover layers for simulation set up.")),
    defineParameter("deciduousCoverDiscount", "numeric", 0.8418911, NA, NA,
                    desc = paste("This was estimated with data from NWT on March 18, 2020 and may or may not be universal.",
                                 "Should be the same as the one used when preparing 'cohortData' in the simulation set up.")),
    defineParameter("LCChangeYr", "integer", c(2001:2011), 1985, 2015,
                    desc = paste("An integer or vector of integers of the validation period years, defining which",
                                 "years of land-cover changes (i.e. disturbances) should be excluded.",
                                 "Only used if rstLCChangeYr is not NULL.",
                                 "See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information.")),
    defineParameter("minCoverThreshold", "numeric", 5, 0, 100,
                    desc = paste("Cover that is equal to or below this number will be omitted from the dataset",
                                 "Should be the same as the one used when preparing 'cohortData' in the simulation set up.")),
    defineParameter("obsDeltaAgeB", "logical", TRUE, NA, NA,
                    desc = paste("When TRUE, the observed changes in biomass and age (deltaB, deltaAge) between",
                                 "the two validation years will be plotted as maps and scatterplots")),
    defineParameter("pixelGroupBiomassClass", "numeric", 100, NA, NA,
                    desc = paste("When assigning pixelGroup membership, this defines the resolution of biomass that will be",
                                 "considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same",
                                 "Should be the same as the one used when preparing 'cohortData' in the simulation set up.")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    desc =  "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("validationReps", "integer", 1:10, NA, NA,
                    desc = paste("The simulation repetitions for the validation. Defaults to 1:10. Set to NA if not using repetitions",
                    "(i.e. only one run)")),
    defineParameter("validationYears", "integer", c(2001, 2011), NA, NA,
                    desc = "The simulation years for the validation. Defaults to 2001 and 2011. Must select two years"),
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA,
                    desc = paste("Vector of length = 1, describing the simulation time at",
                                 "which the first plot event should occur. Set to NA to turn plotting off.")),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, desc = "This describes the simulation time at which the first save event should occur"),
    defineParameter(".savePlots", "logical", TRUE, NA, NA, desc = "Whether plots should be saved in file.path(outputPath(sim), 'Figs')"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, desc = "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default")
  ),
  inputObjects = bind_rows(
    expectsInput("allCohortData", "data.table",
                 desc = paste("All 'cohortData' tables saved during the simulation, particularly for the validation years.",
                              "If not supplied, the module will attempt to retrieve them using the 'simulationOutputs' table")),
    expectsInput("biomassMap", "RasterLayer",
                 desc = paste("total biomass raster layer in study area, filtered for pixels covered by cohortData.",
                              "Only used to calculate total no. of pixels being simulated",
                              "If not supplied, will default to to the Canadian Forestry",
                              "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                              "from 2001. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata")),
    expectsInput("firePerimeters", "sf",
                 desc = paste("A map of fire perimeters in the study area that can be used to exclude pixels",
                              "that have been burnt during the validation period. Defaults to the Canadian",
                              "Wildland Fire Information System 1986-2018 National Burned Area Composite,",
                              "subset to fires between 2001 and 2011 (inclusively)."),
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nbac/nbac_1986_to_2019_20200921.zip"),
    expectsInput("pixelGroupMapStk", "RasterStack",
                 desc = paste("A stack of pixelGroupMaps saved during the simulation, particularly for the validation years.",
                              "If not supplied, the module will attempt to make it using the 'simulationOutputs' table")),
    expectsInput("rawBiomassMapStart", "RasterLayer",
                 desc = paste("observed total biomass raster layer in study area at the first year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived",
                              "total aboveground biomass map from 2001 (in tonnes/ha). If necessary, biomass values",
                              "are rescaled to match changes in resolution.",
                              "See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata."),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")),
    expectsInput("rawBiomassMapEnd", "RasterLayer",
                 desc = paste("observed total biomass raster layer in study area at the last year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived",
                              "total aboveground biomass map from 2011 (in tonnes/ha). If necessary, biomass values",
                              "are rescaled to match changes in resolution.",
                              "See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
                                    "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = paste("A raster of the studyArea in the same resolution and projection as rawBiomassMapStart.",
                              "This is the scale used for all *outputs* for use in the simulation.")),
    expectsInput("rstLCChange", "RasterLayer",
                 desc = paste("A mask-type map of land cover changes in the study area that can be used to exclude pixels",
                              "that have been disturbed during the validation period. Defaults to Canada's forest",
                              "change map between 1985-2011 (CFS), filtered for years 2001-2011 (inclusively)",
                              "and all disturbances collapsed (map only has values of 1 and NA). See parameter LCChangeYr",
                              "to change the period of disturbances, and",
                              "https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information."),
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/C2C_change_year_1985_2011.zip"),
    expectsInput("rstLCChangeYr", "RasterLayer",
                 desc = paste("An OPTIONAL map of land cover change years in the study area used to exclude pixels that have",
                              "been disturbed during the validation period. Defaults to Canada's forest",
                              "change national map between 1985-2011 (CFS). By default disturbances are subset to",
                              " to years 2001-2011 (inclusively; see parameter LCChangeYr).",
                              "See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information."),
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/C2C_change_year_1985_2011.zip"),
    expectsInput("simulationOutputs", "data.table",
                 desc = paste("An OPTIONAL table listing simulation outputs (as passed to 'spades()', or 'experiment",
                              "that will be used to make 'allCohortData', 'pixelGroupMapStk',",
                              "if these are not provided.")),
    expectsInput("speciesLayersStart", "RasterStack",
                 desc = paste("observed cover percentage raster layers by species in Canada species map,",
                              "at the first year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived species cover maps from 2001, using a cover threshold of 10% -",
                              "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/")),
    expectsInput("speciesLayersEnd", "RasterStack",
                 desc = paste("observed cover percentage raster layers by species in Canada species map used for validation",
                              "at the last year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived species cover maps from 2011 -",
                              "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/")),
    expectsInput("sppColorVect", "character",
                 desc = "named character vector of hex colour codes corresponding to each species"),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA."),
    expectsInput("standAgeMapStart", "RasterLayer",
                 desc =  paste("observed stand age map in study area, at the first year of the validation period",
                               "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                               "kNN-derived biomass map from 2001 -",
                               "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("standAgeMapEnd", "RasterLayer",
                 desc = paste("observed stand age raster layer in study area, at the last year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived stand age",
                              "map from 2011. See",
                              "https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
                                    "NFI_MODIS250m_2011_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the study area.",
                              "Defaults to  an area in Southwestern Alberta, Canada.")),
  ),
  outputObjects = bind_rows(
    createsOutput("rstDisturbedPix", "RasterLayer",
                  desc = paste("Raster of pixel IDs (as a mask) that have been disturbed by fire or suffered land-cover",
                               "changes during the validation period. These pixels are excluded form the validation.")),
    createsOutput("rawBiomassMapStart", "RasterLayer",
                  desc = paste("observed total biomass raster layer in study area at the first year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("rawBiomassMapEnd", "RasterLayer",
                  desc = paste("observed total biomass raster layer in study area at the last year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("speciesLayersStart", "RasterStack",
                  desc = paste("observed cover percentage raster layers by species in Canada species map,",
                               "at the first year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("speciesLayersEnd", "RasterStack",
                  desc = paste("observed cover percentage raster layers by species in Canada species map,",
                               "at the last year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("standAgeMapStart", "RasterLayer",
                  desc =  paste("observed stand age map in study area, at the first year of the validation period",
                                "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("standAgeMapEnd", "RasterLayer",
                  desc = paste("observed stand age map in study area, at the last year of the validation period",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("standCohortData", "data.table",
                  desc = paste(""))
  )
))

doEvent.Biomass_validationKNN = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)

      if (!is.na(P(sim)$.plotInitialTime)) {
        dev(height = 7, width = 12)   ## don't overwrite other plots, open new window
        mod$statsWindow <- dev.cur()
        if (P(sim)$obsDeltaAgeB) {
          mod$mapWindow <- mod$statsWindow + 1 ## this window should to be made first if need be
          mod$landscapeWindow <- mod$mapWindow + 1
        } else
          mod$landscapeWindow <- mod$statsWindow + 1

        mod$standWindow <- mod$landscapeWindow + 1
      }

      if (P(sim)$obsDeltaAgeB) {
        sim <- scheduleEvent(sim, eventType = "obsDeltaMaps", eventTime = times(sim)$start,
                             eventPriority = 1, moduleName = currentModule(sim))
      }

      sim <- scheduleEvent(sim, eventType = "landscapeWidePlots", eventTime = times(sim)$start,
                           eventPriority = 2, moduleName = currentModule(sim))
      sim <- scheduleEvent(sim, eventType = "standLevelPlots", eventTime = times(sim)$start,
                           eventPriority = 3, moduleName = currentModule(sim))
      sim <- scheduleEvent(sim, eventType = "deltaBComparisons", eventTime = times(sim)$start,
                           eventPriority = 3, moduleName = currentModule(sim))
    },
    obsDeltaMaps = {
      sim <- obsrvdDeltaMapsEvent(sim)
    },
    landscapeWidePlots = {
      sim <- landscapeWidePlotsEvent(sim)
    },
    standLevelPlots = {
      sim <- standLevelPlotsEvent(sim)
    },
    deltaBComparisons = {
      sim <- deltaBComparisonsEvent(sim)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
Init <- function(sim) {
  ## MAKE PLOT SAVING DIRECTORY IF NEEDED
  if (P(sim)$.savePlots) {
    plotPath <- file.path(outputPath(sim), "Figs")
    if (!dir.exists(plotPath))
      dir.create(plotPath, recursive = TRUE)
    mod$plotPath <- plotPath
  }

  ## make internal module reps variable
  mod$validationReps <- if (!any(is.na(P(sim)$validationReps)))
    P(sim)$validationReps else
      1L

  ## CHECK RASTER LAYERS AGAINST RASTERTOMATCH -----------------------------------
  if (!compareRaster(sim$biomassMap, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$biomassMap <- postProcess(sim$biomassMap, rasterToMatch = sim$rasterToMatch)
  }

  if (!compareRaster(sim$rawBiomassMapStart, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$rawBiomassMapStart <- postProcess(sim$rawBiomassMapStart, rasterToMatch = sim$rasterToMatch)
  }
  if (!compareRaster(sim$rawBiomassMapEnd, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$rawBiomassMapEnd <- postProcess(sim$rawBiomassMapEnd, rasterToMatch = sim$rasterToMatch)
  }

  if (!compareRaster(sim$speciesLayersStart, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$speciesLayersStart <- postProcess(sim$speciesLayersStart, rasterToMatch = sim$rasterToMatch)
  }

  if (!compareRaster(sim$speciesLayersEnd, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$speciesLayersEnd <- postProcess(sim$speciesLayersEnd, rasterToMatch = sim$rasterToMatch)
  }

  if (!compareRaster(sim$standAgeMapStart, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$standAgeMapStart <- postProcess(sim$standAgeMapStart, rasterToMatch = sim$rasterToMatch)
  }

  if (!compareRaster(sim$standAgeMapEnd, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$standAgeMapEnd <- postProcess(sim$standAgeMapEnd, rasterToMatch = sim$rasterToMatch)
  }

  if (!compareRaster(sim$rstLCChange, sim$rasterToMatch, stopiffalse = FALSE)) {
    sim$rstLCChange <- postProcess(sim$rstLCChange, rasterToMatch = sim$rasterToMatch)
  }

  ## EXCLUDE DISTURBED PIXELS FROM VALIDATION  -----------------------------------
  ## make a template raster with IDs
  rasterToMatchIDs <- sim$rasterToMatch
  rasterToMatchIDs <- setValues(rasterToMatchIDs, values = 1:ncell(rasterToMatchIDs))

  ## get pixels inside fire perimeters
  inFireIDs <- getValues(mask(rasterToMatchIDs, sim$firePerimeters))
  inFireIDs <- inFireIDs[!is.na(inFireIDs)] ## faster than na.omit

  ## get pixels inside LCC pixels
  inLCChangeIDs <- rasterToMatchIDs[!is.na(sim$rstLCChange)]

  ## make vector of pixels that are both in fire perimeters and LCChange
  ## export to sim
  sim$disturbedIDs <- union(inFireIDs, inLCChangeIDs)

  sim$rstDisturbedPix <- sim$rasterToMatch
  sim$rstDisturbedPix <- setValues(sim$rstDisturbedPix, values = NA)
  sim$rstDisturbedPix[sim$disturbedIDs] <- 1

  ## exclude these pixels from validation layers
  sim$speciesLayersEnd[sim$disturbedIDs] <- NA
  sim$rawBiomassMapEnd[sim$disturbedIDs] <- NA
  sim$standAgeMapEnd[sim$disturbedIDs] <- NA

  sim$speciesLayersStart[sim$disturbedIDs] <- NA
  sim$rawBiomassMapStart[sim$disturbedIDs] <- NA
  sim$standAgeMapStart[sim$disturbedIDs] <- NA

  ## PREPARE OBSERVED DATA (VALIDATION) PIXEL TABLES   -----------------------------------
  ## need to reproduce some steps in Biomass_borealDataPrep
  ## FIRST VALIDATION
  message(blue("Making pixelCohortData-equivalent tables with observed data: first year"))
  pixelTable <- makePixelTable(speciesLayers = sim$speciesLayersStart,
                               biomassMap = sim$rawBiomassMapStart,
                               standAgeMap = sim$standAgeMapStart,
                               rasterToMatch = sim$rasterToMatch)
  set(pixelTable, NULL, c("initialEcoregionCode", "rasterToMatch"), NULL)
  pixelTable <- unique(pixelTable)

  validationDataStart <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                                   minCoverThreshold = P(sim)$minCoverThreshold)
  validationDataStart <- partitionBiomass(x = P(sim)$deciduousCoverDiscount, validationDataStart)
  set(validationDataStart, NULL, "B",
      asInteger(validationDataStart$B/P(sim)$pixelGroupBiomassClass) * P(sim)$pixelGroupBiomassClass)
  set(validationDataStart, NULL, c("decid", "cover2", "logAge"), NULL)
  set(validationDataStart, NULL, "cover", asInteger(validationDataStart$cover))

  ## remove pixels with missing data
  validationDataStart <- validationDataStart[!is.na(B)]

  ## due to adjustment of B (to nearest 100), totalBiomass has to be recalculated
  validationDataStart[, totalBiomass := sum(B), by = pixelIndex]

  ## LAST VALIDATION YEAR
  message(blue("Making pixelCohortData-equivalent tables with observed data: last year"))
  pixelTable <- makePixelTable(speciesLayers = sim$speciesLayersEnd,
                               biomassMap = sim$rawBiomassMapEnd,
                               standAgeMap = sim$standAgeMapEnd,
                               rasterToMatch = sim$rasterToMatch)
  set(pixelTable, NULL, c("initialEcoregionCode", "rasterToMatch"), NULL)
  pixelTable <- unique(pixelTable)

  validationDataEnd <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                                 minCoverThreshold = P(sim)$minCoverThreshold)
  validationDataEnd <- partitionBiomass(x = P(sim)$deciduousCoverDiscount, validationDataEnd)
  set(validationDataEnd, NULL, "B",
      asInteger(validationDataEnd$B/P(sim)$pixelGroupBiomassClass) * P(sim)$pixelGroupBiomassClass)
  set(validationDataEnd, NULL, c("decid", "cover2", "logAge"), NULL)
  set(validationDataEnd, NULL, "cover", asInteger(validationDataEnd$cover))

  ## remove pixels with missing data
  validationDataEnd <- validationDataEnd[!is.na(B)]

  ## due to adjustment of B (to nearest 100), totalBiomass has to be recalculated
  validationDataEnd[, totalBiomass := sum(B), by = pixelIndex]

  ## PREPARE SIMULATED DATA (COHORTDATA) PIXEL TABLES   -----------------------------------
  ## iterate through pixelGroupMaps and subset allCohortData per year/rep to expand the tables
  ## no need to save to sim$
  message(blue("Making pixelCohortData tables of simulated outputs"))

  ## check that starting conditions are the same
  assertRepsAllCohortData(allCohortData = sim$allCohortData, reps = mod$validationReps,
                          years = P(sim)$validationYears)
  allPixelCohortData <- rbindlist(fill = TRUE, use.names = TRUE,
                                  l = lapply(unstack(sim$pixelGroupMapStk), FUN = function(pixelGroupMap, allCohortData) {
                                    yr <- as.numeric(sub("year([[:digit:]]*)_rep.*", "\\1", names(pixelGroupMap)))
                                    rp <- as.numeric(sub(".*_rep", "", names(pixelGroupMap)))
                                    pixelCohortData <- allCohortData[year == yr & rep == rp]
                                    pixelCohortData <- addPixels2CohortData(pixelCohortData, pixelGroupMap)
                                    pixelCohortData
                                  }, allCohortData = sim$allCohortData))

  ## summarize allPixelCohortData to stand totalB per species and
  ## and biomass-averaged stand age
  standCohortData <- allPixelCohortData[, .(B, sum(B), age),
                                        by = .(rep, year, pixelIndex, speciesCode)]
  standCohortData <- standCohortData[, standAge := sum(age * B, na.rm = TRUE) / sum(B, na.rm = TRUE),
                                     by = .(rep, year, pixelIndex)]
  ## drop unnecessary columns and remove separate cohorts
  standCohortData[, B := V2] ## overwrite
  standCohortData[, `:=`(V2 = NULL, age = NULL)]
  standCohortData <- unique(standCohortData)

  ## MERGE OBSERVED AND SIMULATED DATA TABLES  -----------------------------------
  ## add observed data to standCohortData
  ## note that some pixelIndex X spp combinations are lacking because the observed data
  ## has spp in some pixels that are not found in the simulation data,
  ## and vice-versa. To make sure that the observed pixel X spp combinations are added to each
  ## rep/year the tables need to be extended - otherwise sometimes the observed data
  ## is only joined to some reps/years, making observed averages "vary" across reps/years
  combinationsStart <- as.data.table(expand.grid(list(speciesCode = unique(standCohortData$speciesCode),
                                                      pixelIndex = unique(c(validationDataStart$pixelIndex, standCohortData$pixelIndex)),
                                                      rep = mod$validationReps,
                                                      year = P(sim)$validationYears[1])))

  combinationsEnd <- as.data.table(expand.grid(list(speciesCode = unique(standCohortData$speciesCode),
                                                    pixelIndex = unique(c(validationDataEnd$pixelIndex, standCohortData$pixelIndex)),
                                                    rep = mod$validationReps,
                                                    year = P(sim)$validationYears[2])))
  validationDataStart <- validationDataStart[combinationsStart,
                                             on = c("pixelIndex", "speciesCode")]
  validationDataEnd <- validationDataEnd[combinationsEnd,
                                         on = c("pixelIndex", "speciesCode")]
  ## merge and export to sim
  sim$validationData <- rbindlist(list(validationDataStart, validationDataEnd),
                                  use.names = TRUE)

  ## exclude pixels that are not simulated
  sim$validationData <- sim$validationData[pixelIndex %in% standCohortData$pixelIndex]

  ## convert NAs to 0s
  cols <- c("cover", "age", "B", "totalBiomass")
  sim$validationData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]

  ## change names before joining.
  setnames(sim$validationData,
           old = c("cover", "age", "B", "totalBiomass"),
           new = c("coverObsrvd", "standAgeObsrvd", "BObsrvd", "standBObsrvd"))
  ## reorder columns
  cols <- c("rep", "year", "pixelIndex", "speciesCode", "coverObsrvd",
            "standAgeObsrvd", "BObsrvd", "standBObsrvd")
  sim$validationData <- sim$validationData[, ..cols]

  ## merge and keep all combos
  standCohortData <- standCohortData[sim$validationData,
                                     on = c("rep", "year", "pixelIndex", "speciesCode")]

  ## remove disturbed pixels
  standCohortData <- standCohortData[!pixelIndex %in% sim$disturbedIDs]

  ## convert NAs to 0s
  cols <- c("standAge", "B")
  standCohortData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]

  ## calculate simulated standAge, standB and relative B
  standCohortData[, standB := asInteger(sum(B)),
                  .(rep, year, pixelIndex)]
  standCohortData[, relativeAbund := B/standB]
  standCohortData[standB == 0, relativeAbund := 0]

  ## replace inserted standB/AgeObsrvd 0s (coming from merge) with actual stand biomass/age
  standCohortData[, `:=`(standBObsrvd = max(standBObsrvd),
                         standAgeObsrvd = max(standAgeObsrvd)),
                  by = .(rep, year, pixelIndex)]

  ## calculate observed relative B
  standCohortData[, relativeAbundObsrvd := BObsrvd/standBObsrvd]
  standCohortData[standBObsrvd == 0, relativeAbundObsrvd := 0]

  ## classify pixels by dominant species (i,e, veg type) - will need to be corrected for compeeting dominants
  ## these match with model outputs, I checked
  ## note2: mixed pixels get a "mixed" type
  standCohortData[standB > 0, vegType := speciesCode[which.max(relativeAbund)],
                  by = .(year, rep, pixelIndex)]
  standCohortData[standBObsrvd > 0, vegTypeObsrvd := speciesCode[which.max(relativeAbundObsrvd)],
                  by = .(year, rep, pixelIndex)]


  ## get number of dominant species -- note that stands with 0 B will appear as mixed initially, but are fixed below
  tempDT <- standCohortData[, list(noDoms = sum(relativeAbund == max(relativeAbund)),
                                   noDomsObsrvd = sum(relativeAbundObsrvd == max(relativeAbundObsrvd))),
                            by = .(year, rep, pixelIndex)]

  standCohortData <- tempDT[standCohortData, on = .(year, rep, pixelIndex)]
  standCohortData[standB == 0, noDoms := 0]
  standCohortData[standBObsrvd == 0, noDomsObsrvd := 0]

  standCohortData[noDoms > 1, vegType := "Mixed"]
  standCohortData[noDomsObsrvd > 1, vegTypeObsrvd := "Mixed"]

  standCohortData[is.na(vegType), vegType := "No veg."]
  standCohortData[is.na(vegTypeObsrvd), vegTypeObsrvd := "No veg."]

  ## calculate some landscape metrics
  standCohortData[, `:=`(landscapeB = sum(B),
                         landscapeBObsrvd = sum(BObsrvd)),
                  by = .(rep, year)]

  ## reorder column names
  cols <- c(grep("Obsrvd", names(standCohortData), value = TRUE, invert = TRUE),
            grep("Obsrvd", names(standCohortData), value = TRUE))
  standCohortData <- standCohortData[, ..cols]

  ## assert and export to sim -- vegType cols can have NAs
  assertStandCohortData(standCohortData)
  sim$standCohortData <- standCohortData

  ## make labels for plots
  sim$speciesLabels <- equivalentName(unique(sim$standCohortData$speciesCode), sim$sppEquiv,
                                      column = "EN_generic_short")
  names(sim$speciesLabels) <- unique(sim$standCohortData$speciesCode)

  ## EXCLUDE PIXELS WHERE OBSERVED STAND B OR STAND AGE DECREASED -----------------
  ## Only keep pixels where stand age AND stand B increased, or remained the same,
  ## as we cannot account for disturbances that may not have been captured from sat data
  ## or measurement errors that yielded too high B in 2001. It is unlikely that in 10 years we have a large proportion
  ## of stands seeing reduced stand age and B due to death from long-age.
  ## THIS HAS BEEN DESACTIVATED FOR NOW.

  if (FALSE) {
    year1 <- P(sim)$validationYears[1]
    year2 <- P(sim)$validationYears[2]
    tempDT <- unique(sim$standCohortData[, .(year, rep, pixelIndex, standAgeObsrvd, standBObsrvd)])
    tempDT <- tempDT[, .(standDeltaAgeObsrvd = standAgeObsrvd[year == year2] - standAgeObsrvd[year == year1],
                         standDeltaBObsrvd = standBObsrvd[year == year2] - standBObsrvd[year == year1]),
                     by = .(rep, pixelIndex)]
    pixToKeep <- unique(tempDT[standDeltaBObsrvd >= 0 & standDeltaAgeObsrvd >= 0,
                               pixelIndex])
  } else
    pixToKeep <- unique(sim$standCohortData$pixelIndex)

  ## return some statistics about excluded pixels
  pixToRm <- unique(c(setdiff(unique(sim$standCohortData$pixelIndex), pixToKeep),
                      sim$disturbedIDs))
  excludedPixStats <- data.table(noPixels = length(pixToRm),
                                 landscapePrc = round(length(pixToRm)/
                                                        sum(!is.na(getValues(sim$biomassMap))),
                                                      2) * 100)
  message(blue("Pixels disturbed during the validation period, and pixels that showed decreases in age and biomass
               will be excluded from validation, representing a loss of:", excludedPixStats$noPixels, "pixels or",
               excludedPixStats$landscapePrc, "% of the initial simulated landscape."))

  ## keep aforementioned pixels only
  sim$standCohortData <- sim$standCohortData[pixelIndex %in% pixToKeep]

  ## clean up and free memory
  rm(pixelTable, standCohortData, combinationsStart, combinationsEnd, validationDataStart, validationDataEnd)
  .gc()

  return(invisible(sim))
}

obsrvdDeltaMapsEvent <- function(sim) {
  ## MAPS OF OBSERVED CHANGES IN STAND B AND AGE - RAW DATA -------------------
  standDeltaBObsrvdRas <- (sim$rawBiomassMapEnd - sim$rawBiomassMapStart) * 100 ## to rescale to to/ha
  pixToNA <- setdiff(1:ncell(standDeltaBObsrvdRas), sim$standCohortData$pixelIndex)
  standDeltaBObsrvdRas[pixToNA] <- NA

  standDeltaAgeObsrvdRas <- sim$standAgeMapEnd - sim$standAgeMapStart
  standDeltaAgeObsrvdRas[pixToNA] <- NA

  ## what is the relationship between the two?
  standDeltaObsrvdData <- na.omit(data.table(pixelIndex = 1:ncell(standDeltaBObsrvdRas),
                                             standDeltaBObsrvd = getValues(standDeltaBObsrvdRas),
                                             standDeltaAgeObsrvd = getValues(standDeltaAgeObsrvdRas)))

  plot1 <- ggplot(standDeltaObsrvdData,
                  aes(x = standDeltaAgeObsrvd, y = standDeltaBObsrvd)) +
    geom_point() +
    stat_smooth(method = "lm") +
    theme_pubr(base_size = 12, margin = FALSE) +
    labs(y = expression(paste("observed ", Delta, "B"))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())

  plot2 <- ggplot(standDeltaObsrvdData,
                  aes(y = standDeltaBObsrvd)) +
    geom_boxplot() +
    theme_pubr(base_size = 12, margin = FALSE) +
    theme(axis.title.y = element_blank(), axis.text = element_blank())

  plot3 <- ggplot(standDeltaObsrvdData,
                  aes(y = standDeltaAgeObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste("observed ", Delta, "age"))) +
    theme_pubr(base_size = 12, margin = FALSE) +
    theme(axis.text.y = element_blank())

  plot4 <- ggarrange(plot1, plot2,  plot3, nrow = 2, ncol = 2, align = "hv",
                     widths = c(1, 0.5), heights = c(1, 0.5))

  ## delta biomass for the "supposed" age increment - all over the place
  year1 <- P(sim)$validationYears[1]
  year2 <- P(sim)$validationYears[2]
  yearGap <- year2 - year1
  plot5 <- ggplot(standDeltaObsrvdData[standDeltaAgeObsrvd == yearGap],
                  aes(y = standDeltaBObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste("observed ", Delta, "B")),
         title = bquote(atop("observed" ~ Delta ~ "B", "in pixels that aged" ~ .(yearGap)))) +
    theme_pubr(base_size = 12, margin = FALSE)

  ## MAPS OF OBSERVED CHANGES IN STAND B AND AGE - AFTER ADJUSTMENTS -------
  ## by adjustments we mean, the data cleanup that we replicate from Biomass_borealDataPrep
  plotData <- unique(sim$standCohortData[, .(year, pixelIndex, standBObsrvd, standAgeObsrvd)])
  plotData <- plotData[, list(standDeltaBObsrvd = unique(standBObsrvd[which(year == year2)]) - unique(standBObsrvd[which(year == year1)]),
                              standDeltaAgeObsrvd = unique(standAgeObsrvd[which(year == year2)]) - unique(standAgeObsrvd[which(year == year1)])),
                       , by = pixelIndex]

  if (any(duplicated(plotData$pixelIndex)))
    stop("There should not be duplicated pixels in observed data")

  standDeltaBObsrvdAdj <- setValues(sim$rasterToMatch, NA)
  standDeltaBObsrvdAdj[plotData[, pixelIndex]] <- plotData[, standDeltaBObsrvd]

  standDeltaAgeObsrvdAdj <- setValues(sim$rasterToMatch, NA)
  standDeltaAgeObsrvdAdj[plotData[, pixelIndex]] <- plotData[, standDeltaAgeObsrvd]


  plot6 <- ggplot(plotData,
                  aes(x = standDeltaAgeObsrvd, y = standDeltaBObsrvd)) +
    geom_point() +
    stat_smooth(method = "lm") +
    theme_pubr(base_size = 12, margin = FALSE) +
    labs(y = expression(paste("observed ", Delta, "B", " - adjusted"))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())

  plot7 <- ggplot(standDeltaObsrvdData,
                  aes(y = standDeltaBObsrvd)) +
    geom_boxplot() +
    theme_pubr(base_size = 12, margin = FALSE) +
    theme(axis.title.y = element_blank(), axis.text = element_blank())

  plot8 <- ggplot(standDeltaObsrvdData,
                  aes(y = standDeltaAgeObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste("observed ", Delta, "age", " - adjusted"))) +
    theme_pubr(base_size = 12, margin = FALSE) +
    theme(axis.text.y = element_blank())

  plot9 <- ggarrange(plot6, plot7,  plot8, nrow = 2, ncol = 2, align = "hv",
                     widths = c(1, 0.5), heights = c(1, 0.5))

  ## delta biomass for the "supposed" age increment - all over the place
  plot10 <- ggplot(plotData[standDeltaAgeObsrvd == yearGap],
                   aes(y = standDeltaBObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste("observed ", Delta, "B", " - adjusted")),
         title = bquote(atop("observed" ~ Delta ~ "B - adjusted", "in pixels that aged" ~ .(yearGap)))) +
    theme_pubr(base_size = 12, margin = FALSE)

  if (!is.na(P(sim)$.plotInitialTime)) {
    dev(mod$mapWindow)
    clearPlot()
    Plot(standDeltaBObsrvdRas, standDeltaAgeObsrvdRas, new = TRUE,
         title = c("stand delta-B", "stand delta-Age"))
    Plot(standDeltaBObsrvdAdj, standDeltaAgeObsrvdAdj, new = TRUE,
         title = c("stand delta-B - adjusted", "stand delta-Age - adjusted"))

    dev(mod$statsWindow)
    clearPlot()
    plotDeltaStats <- ggarrange(plot4, plot5, plot9, plot10)
    Plot(plotDeltaStats, new = TRUE, title = "") ## title = FALSE not working
  }

  if (P(sim)$.savePlots) {
    ggsave(filename = file.path(mod$plotPath, "observedDeltaBDeltaAge_lm.png"),
           plot = plot4, width = 7, height = 5, units = "in")

    ggsave(filename = file.path(mod$plotPath, "observedDeltaB_yearGap.png"),
           plot = plot5, width = 5, height = 4, units = "in")

    ggsave(filename = file.path(mod$plotPath, "observedDeltaBDeltaAge_lmADJ.png"),
           plot = plot9, width = 7, height = 5, units = "in")

    ggsave(filename = file.path(mod$plotPath, "observedDeltaB_yearGapADJ.png"),
           plot = plot10, width = 5, height = 4, units = "in")

    if (!is.na(P(sim)$.plotInitialTime)) {
      dev.set(mod$mapWindow)
      dev.copy(pdf, file = file.path(mod$plotPath, 'deltaB_Age_Maps.pdf'))
      dev.off()
    } else {
      pdf(file = file.path(mod$plotPath, 'deltaB_Age_Maps.pdf'))
      layout(matrix(1:4, ncol=2))
      raster::plot(standDeltaBObsrvdRas, main = "stand delta-B")
      raster::plot(standDeltaAgeObsrvdRas, main = "stand delta-Age")
      raster::plot(standDeltaBObsrvdAdj, main = "stand delta-B - adjusted")
      raster::plot(standDeltaAgeObsrvdAdj, main = "stand delta-Age - adjusted")
      dev.off()
    }

  }

  return(invisible(sim))
}

landscapeWidePlotsEvent <- function(sim) {
  assertStandCohortData(sim$standCohortData)

  plotData <- sim$standCohortData
  plotData <- plotData[, list(landRelativeAbund = sum(B)/unique(landscapeB),
                              landRelativeAbundObsrvd = sum(BObsrvd)/unique(landscapeBObsrvd)),
                       by = .(rep, year, speciesCode)]

  plot1 <- ggplot(data = plotData,
                  aes(x = speciesCode, y = landRelativeAbund)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    stat_summary(aes(y = landRelativeAbundObsrvd, colour = "observed"),
                 fun = "mean", geom = "point", size = 2) +
    scale_x_discrete(labels = sim$speciesLabels, drop = FALSE) +
    scale_color_manual(values = c("observed" = "red3")) +
    theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
    theme(legend.position = "right") +
    facet_wrap(~ year) +
    labs(title = "Species relative abundances",
         x = "", y = expression(over("species B", "total B")),
         colour = "")

  ## no. pixels with a species
  plotData <- sim$standCohortData
  plotData <- plotData[,  list(count = sum(B > 0),
                               countObsrvd = sum(BObsrvd > 0)),
                       by = .(rep, year, speciesCode)]
  plotData <- melt.data.table(plotData, measure.vars = c("count", "countObsrvd"),
                              variable.name = "dataType", value.name = "count")

  plot2 <- ggplot(data = plotData[dataType == "count"],
                  aes(x = speciesCode, y = count)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    stat_summary(data = plotData[dataType == "countObsrvd"],
                 aes(x = speciesCode, y = count, colour = "observed"),
                 fun = "mean", geom = "point", size = 2) +
    scale_x_discrete(labels = sim$speciesLabels, drop = FALSE) +
    scale_color_manual(values = c("observed" = "red3")) +
    theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
    theme(legend.position = "right") +
    facet_wrap(~ year) +
    labs(title = "Species presences", x = "", y = "no. pixels",
         colour = "", fill = "")

  ## no. pixels with a certain dominant species
  ## note: don't use melt, because dominant spp differ between valid and simul data.
  ## simulated and observed differ in no. of pixels in year 1...
  ## this is because B is adjusted using a statistical model
  plotData1 <- unique(sim$standCohortData[, .(year, rep, pixelIndex, vegType)])
  plotData1[, dataType := "simulated"]
  plotData2 <- unique(sim$standCohortData[, .(year, rep, pixelIndex, vegTypeObsrvd)])
  plotData2[, dataType := "observed"]
  setnames(plotData2, "vegTypeObsrvd", "vegType")

  plotData <- rbind(plotData1, plotData2)
  rm(plotData1, plotData2)


  plotData <- plotData[, list(count = length(pixelIndex)),
                        by = .(rep, year, dataType, vegType)]
  plot3 <- ggplot(data = plotData[dataType == "simulated"],
                  aes(x = vegType, y = count)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    geom_point(data = plotData[dataType == "observed"],
               aes(x = vegType, y = count, colour = "observed"), size = 2) +
    scale_x_discrete(labels = sim$speciesLabels, drop = FALSE) +
    scale_color_manual(values = c("observed" = "red3")) +
    theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
    theme(legend.position = "right") +
    facet_wrap(~ year) +
    labs(title = "Dominant species' presences",
         x = "", y = "no. pixels", fill = "", colour = "")

  maxPixels <- sum(!is.na(getValues(sim$biomassMap)))
  plotLandscapeComp <- ggarrange(plot1 + scale_y_continuous(limits = c(0,1)),
                                 plot2 + scale_y_continuous(limits = c(0, maxPixels)),
                                 plot3 + scale_y_continuous(limits = c(0, maxPixels)),
                                 common.legend = TRUE, legend = "bottom",
                                 nrow = 2, ncol = 2)

  if (!is.na(P(sim)$.plotInitialTime)) {
    dev(mod$landscapeWindow)
    clearPlot()
    Plot(plotLandscapeComp, title = "Landscape-averaged comparisons", new = TRUE)
  }

  if (P(sim)$.savePlots) {
    plotLandscapeComp2 <- annotate_figure(plotLandscapeComp,
                                          top = text_grob("Landscape-averaged comparisons", size = 16))
    ggsave(filename = file.path(mod$plotPath, "LandscapeComparisons_relB_PresAbs.png"),
           plot = plotLandscapeComp2, width = 12, height = 7, units = "in")
  }

  return(invisible(sim))
}

standLevelPlotsEvent <- function(sim) {
  ## STAND/PIXEL-LEVEL COMPARISONS IN A GIVEN YEAR --------------------
  assertStandCohortData(sim$standCohortData)

  ## stand/pixel-level relative abundances per species
  plotData <- sim$standCohortData
  plotData <- plotData[, .(rep, year, pixelIndex, speciesCode,
                           relativeAbund, relativeAbundObsrvd)]
  plotData <- melt.data.table(plotData,
                              id.vars = c("rep", "year", "pixelIndex", "speciesCode"))
  plot1 <- ggplot(data = plotData,
                  aes(x = speciesCode, y = value, fill = variable)) +
    geom_boxplot() +
    scale_x_discrete(labels = sim$speciesLabels, drop = FALSE) +
    scale_fill_discrete(labels = c("relativeAbund" = "simulated",
                                   "relativeAbundObsrvd" = "observed")) +
    theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
    facet_wrap(~ year) +
    labs(title = "Species relative abundances", fill = "",
         x = "", y = expression(over("species B", "stand B")))

  ## stand/pixel-level relative abundance per dominant species
  ## get dominant species (those with maxB) - these match with model outputs, I checked
  tempDT <- sim$standCohortData
  tempDT <- tempDT[, list(noDoms = sum(relativeAbund == max(relativeAbund)),
                          noDomsObsrvd = sum(relativeAbundObsrvd == max(relativeAbundObsrvd))),
                   by = .(year, rep, pixelIndex)]

  plotData <- sim$standCohortData
  plotData <- unique(plotData[, list(vegType = speciesCode[which.max(relativeAbund)],
                                     relativeAbund = max(relativeAbund),
                                     vegTypeObsrvd = speciesCode[which.max(relativeAbundObsrvd)],
                                     relativeAbundObsrvd = max(relativeAbundObsrvd)),
                              by = .(year, rep, pixelIndex)])
  plotData <- tempDT[plotData, on = .(year, rep, pixelIndex)]

  plotData1 <- unique(plotData[, .(year, rep, pixelIndex, vegType, relativeAbund, noDoms)])
  plotData1[, dataType := "relativeAbund"]
  plotData2 <- unique(plotData[, .(year, rep, pixelIndex, vegTypeObsrvd, relativeAbundObsrvd, noDomsObsrvd)])
  plotData2[, dataType := "relativeAbundObsrvd"]
  setnames(plotData2, c("vegTypeObsrvd", "relativeAbundObsrvd", "noDomsObsrvd"),
           c("vegType", "relativeAbund", "noDoms"))

  plotData <- rbind(plotData1, plotData2)
  plotData[noDoms > 1, vegType := "Mixed"]
  rm(plotData1, plotData2)

  plot2 <- ggplot(data = plotData,
                  aes(x = vegType, y = relativeAbund, fill = dataType)) +
    geom_boxplot() +
    scale_fill_discrete(labels = c("relativeAbund" = "simulated",
                                   "relativeAbundObsrvd" = "observed")) +
    scale_x_discrete(labels = sim$speciesLabels, drop = FALSE) +
    theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
    facet_wrap(~ year) +
    labs(title = "Dominant species' relative abundances",
         x = "", y = expression(over("species B", "total B")),
         fill = "")

  standCompPlot <- ggarrange(plot1,
                             plot2 + labs(y = " \n "),
                             common.legend = TRUE, legend = "bottom",
                             ncol = 2)

  if (!is.na(P(sim)$.plotInitialTime)) {
    dev(mod$standWindow)
    clearPlot()
    Plot(standCompPlot, title = "Stand-level comparisons", new = TRUE)
  }

  if (P(sim)$.savePlots) {
    standCompPlot2 <- annotate_figure(standCompPlot,
                                      top = text_grob("Stand-level comparisons", size = 16))
    ggsave(filename = file.path(mod$plotPath, "StandComparisons_relB.png"),
           plot = standCompPlot2, width = 12, height = 6, units = "in")
  }

  return(invisible(sim))
}

deltaBComparisonsEvent <- function(sim) {
  ## COMPARISONS OF DELTA B PER PIXEL-------------------
  ## per species
  year1 <- P(sim)$validationYears[1]
  year2 <- P(sim)$validationYears[2]

  plotData <- sim$standCohortData
  plotData <- plotData[, list(deltaB = as.numeric(B[which(year == year2)] - B[which(year == year1)]),
                              deltaBObsrvd = as.numeric(BObsrvd[which(year == year2)] - BObsrvd[which(year == year1)]),
                              standDeltaB = as.numeric(standB[which(year == year2)] - standB[which(year == year1)]),
                              standDeltaBObsrvd = as.numeric(standBObsrvd[which(year == year2)] - standBObsrvd[which(year == year1)])),
                       by = .(rep, pixelIndex, speciesCode)]

  ## melt spp and stand delta separately and rbind
  cols <- c("rep", "pixelIndex", "speciesCode", "deltaB", "deltaBObsrvd")
  plotData1 <- melt.data.table(plotData[, ..cols], measure.vars = c("deltaB", "deltaBObsrvd"),
                               variable.name = "dataType", value.name = "deltaB")
  cols <- c("rep", "pixelIndex", "speciesCode", "standDeltaB", "standDeltaBObsrvd")
  plotData2 <- melt.data.table(plotData[, ..cols], measure.vars = c("standDeltaB", "standDeltaBObsrvd"),
                               variable.name = "dataType", value.name = "deltaB")
  plotData2[dataType == "standDeltaB", dataType := "deltaB"]
  plotData2[dataType == "standDeltaBObsrvd", dataType := "deltaBObsrvd"]
  plotData2 <- unique(plotData2[, speciesCode := "stand"])

  plotData <- rbind(plotData1, plotData2, use.names = TRUE)
  rm(plotData1, plotData2)

  plot1 <-  ggplot(data = plotData[dataType == "deltaB"],
                   aes(x = speciesCode, y = deltaB, group = rep)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    stat_summary(data = plotData[dataType == "deltaBObsrvd"],
                 aes(x = speciesCode, y = deltaB, group = rep),
                 fun = "mean", geom = "point", size = 2, colour = "red3") +
    scale_x_discrete(labels = sim$speciesLabels, drop = FALSE) +
    theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
    labs(title = "Landscape-averaged",
         x = "", y = expression(paste(Delta, "B")))

  plot2 <- ggplot(data = plotData,
                  aes(x = speciesCode, y = deltaB, fill = dataType)) +
    geom_boxplot(aes(alpha = speciesCode == "stand")) +
    scale_x_discrete(labels = c(sim$speciesLabels, "stand" = "Stand")) +
    scale_fill_discrete(labels = c("deltaB" = "simulated",
                                   "deltaBObsrvd" = "observed")) +
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 1.0), guide = FALSE) +
    theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
    labs(title = "Stand-level", fill = "",
         x = "", y = expression(paste(Delta, "B")))

  simObsDeltaBPlot <- ggarrange(plot1, plot2 + labs(y = " \n "),
                                common.legend = TRUE, legend = "bottom",
                                ncol = 2)

  if (!is.na(P(sim)$.plotInitialTime)) {
    dev.set(mod$statsWindow)
    Plot(simObsDeltaBPlot, title = "", new = TRUE)
  }

  if (P(sim)$.savePlots) {
    simObsDeltaBPlot2 <- annotate_figure(simObsDeltaBPlot,
                                         top = text_grob("Stand-level comparisons", size = 16))
    ggsave(filename = file.path(mod$plotPath, "LandscapeStandComparisons_deltaB.png"),
           plot = simObsDeltaBPlot2, width = 10, height = 6, units = "in")
  }
  return(invisible(sim))
}

## INPUT OBJECTS ------------------------------------------
.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

  ## check that there are two validation years
  if (length(P(sim)$validationYears) != 2)
    stop("You must pass TWO 'validationYears' (e.g. 'c(2001L, 2011L)').")


  ## Study area raster and shapefile -------------------------------------------
  if (!suppliedElsewhere("studyArea", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)
  }

  needRTM <- FALSE
  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      needRTM <- TRUE
      message("There is no rasterToMatch supplied; will attempt to use rawBiomassMapStart")
    } else {
      stop("rasterToMatch is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (!suppliedElsewhere("rawBiomassMapStart", sim) || needRTM) {
    rawBiomassMapFilename <- "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
    httr::with_config(config = httr::config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
      #necessary for KNN
      sim$rawBiomassMapStart <- Cache(prepInputs,
                                      url = extractURL("rawBiomassMapStart"),
                                      destinationPath = dPath,
                                      studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMapStart, LCC.. etc
                                      rasterToMatch = if (!needRTM) sim$rasterToMatchLarge else NULL,
                                      maskWithRTM = if (!needRTM) TRUE else FALSE,
                                      useSAcrs = FALSE,     ## never use SA CRS
                                      method = "bilinear",
                                      datatype = "INT2U",
                                      filename2 = .suffix("rawBiomassMapStart.tif", paste0("_", P(sim)$.studyAreaName)),
                                      overwrite = TRUE,
                                      userTags = c(cacheTags, "rawBiomassMapStart"),
                                      omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    })

    ## if using custom raster resolution, need to allocate biomass proportionally to each pixel
    ## if no rawBiomassMapStart/RTM/RTMLarge were suppliedElsewhere, the "original" pixel size respects
    ## whatever resolution comes with the rawBiomassMapStart data
    simPixelSize <- unique(asInteger(res(sim$rasterToMatchLarge)))
    origPixelSize <- 250L # unique(res(sim$rawBiomassMapStart)) ## TODO: figure out a good way to not hardcode this

    if (simPixelSize != origPixelSize) { ## make sure we are comparing integers, else else %!=%
      rescaleFactor <- (origPixelSize / simPixelSize)^2
      sim$rawBiomassMapStart <- sim$rawBiomassMapStart / rescaleFactor
    }
  }

  if (needRTM) {
    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMapStart
    sim$rasterToMatch <- sim$rawBiomassMapStart
    RTMvals <- getValues(sim$rasterToMatch)
    sim$rasterToMatch[!is.na(RTMvals)] <- 1

    sim$rasterToMatch <- Cache(writeOutputs, sim$rasterToMatch,
                               filename2 = file.path(cachePath(sim), "rasters", "rasterToMatch.tif"),
                               datatype = "INT2U", overwrite = TRUE,
                               userTags = c(cacheTags, "rasterToMatch"),
                               omitArgs = c("userTags"))

    ## this is old, and potentially not needed anymore
    if (FALSE) {
      studyArea <- sim$studyArea # temporary copy because it will be overwritten if it is suppliedElsewhere
      message("  Rasterizing the studyArea polygon map")
      if (!is(studyArea, "SpatialPolygonsDataFrame")) {
        dfData <- if (is.null(rownames(studyArea))) {
          polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
          data.frame("field" = as.character(seq_along(length(studyArea))), row.names = polyID)
        } else {
          polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
          data.frame("field" = rownames(studyArea), row.names = polyID)
        }
        studyArea <- SpatialPolygonsDataFrame(studyArea, data = dfData)
      }
      if (!identical(crs(studyArea), crs(sim$rasterToMatch))) {
        studyArea <- spTransform(studyArea, crs(sim$rasterToMatch))
        studyArea <- fixErrors(studyArea)
      }


      #TODO: review whether this is necessary (or will break LandWeb if removed) see Git Issue #22
      # layers provided by David Andison sometimes have LTHRC, sometimes LTHFC ... chose whichever
      LTHxC <- grep("(LTH.+C)", names(studyArea), value = TRUE)
      fieldName <- if (length(LTHxC)) {
        LTHxC
      } else {
        if (length(names(studyArea)) > 1) {
          ## study region may be a simple polygon
          names(studyArea)[1]
        } else NULL
      }

      sim$rasterToMatch <- crop(fasterizeFromSp(studyArea, sim$rasterToMatch, fieldName),
                                studyArea)
      sim$rasterToMatch <- Cache(writeRaster, sim$rasterToMatch,
                                 filename = file.path(dPath, "rasterToMatch.tif"),
                                 datatype = "INT2U", overwrite = TRUE,
                                 userTags = c(cacheTags, "rasterToMatch"),
                                 omitArgs = c("userTags"))
    }
  }

  if (!identical(crs(sim$studyArea), crs(sim$rasterToMatch))) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  ## Land cover change (type and year) -------------------------------------------
  if (!suppliedElsewhere("rstLCChange", sim)) {
    ## get pixels that suffered land-cover changes during the
    ## validation period
    ## change codes are
    ## 0 (or NA?) = no change
    ## 1 = wildfire
    ## 2 = harvest
    ## 5 = low confidence wildfire
    ## 6 = low confidence harvest
    ## change years are from 1985-2010, but values range 85-110

    ## need the year of change map to subset CFSs land-cover change type map.
    LCChangeFilename <- "C2C_change_type_1985_2011.tif"
    LCChangeYrFilename <- "C2C_change_year_1985_2011.tif"
    sim$rstLCChange <- Cache(prepInputs,
                             targetFile = LCChangeFilename,
                             archive = asPath("C2C_change_type_1985_2011.zip"),
                             url = extractURL("rstLCChange"),
                             destinationPath = dPath,
                             studyArea = sim$studyArea,
                             rasterToMatch = sim$rasterToMatch,
                             useSAcrs = FALSE,
                             maskWithRTM = TRUE,
                             method = "ngb",
                             datatype = "INT2U",
                             filename2 = TRUE, overwrite = TRUE,
                             userTags = c("rstLCChange", cacheTags),
                             omitArgs = c("destinationPath", "targetFile", "userTags"))
    ## convert to mask
    sim$rstLCChange[!is.na(sim$rstLCChange)] <- 1

    sim$rstLCChangeYr <- Cache(prepInputs,
                               targetFile = LCChangeYrFilename,
                               archive = asPath("C2C_change_year_1985_2011.zip"),
                               url = extractURL("rstLCChangeYr"),
                               destinationPath = dPath,
                               studyArea = sim$studyArea,
                               rasterToMatch = sim$rasterToMatch,
                               useSAcrs = FALSE,
                               maskWithRTM = TRUE,
                               method = "ngb",
                               datatype = "INT2U",
                               filename2 = TRUE, overwrite = TRUE,
                               userTags = c("rstLCChangeYr", cacheTags),
                               omitArgs = c("destinationPath", "targetFile", "userTags"))

    ## only keep pixels that have been disturbed during the validation period
    ## convert years to the map's format
    yrs <- P(sim)$LCChangeYr - 1900
    pixKeep <- !is.na(getValues(sim$rstLCChange)) &
      getValues(sim$rstLCChangeYr) %in% yrs

    sim$rstLCChange[!pixKeep] <- NA
    sim$rstLCChangeYr[!pixKeep] <- NA
  }

  ## Check that rstLCChange is a mask and matches RTM
  assertRstLCChange(sim$rstLCChange, sim$rasterToMatch)

  ## Fire perimeter data ---------------------------------------------------

  if (!suppliedElsewhere("firePerimeters", sim)) {
    firePerimetersFile <- "nbac_1986_to_2019_20200921.shp"
    sim$firePerimeters <- Cache(prepInputs,
                                targetFile = firePerimetersFile,
                                alsoExtract = "similar",
                                archive = asPath("nbac_1986_to_2019_20200921.zip"),
                                url = extractURL("firePerimeters"),
                                destinationPath = dPath,
                                studyArea = sim$studyArea,
                                rasterToMatch = sim$rasterToMatch,
                                useSAcrs = FALSE,
                                datatype = "INT2U",
                                filename2 = TRUE, overwrite = TRUE,
                                userTags = c("firePerimeters", cacheTags),
                                omitArgs = c("destinationPath", "targetFile", "userTags"))
    ## convert to sf
    sim$firePerimeters <- st_as_sf(sim$firePerimeters)

    ## exclude fire years outside validation period
    sim$firePerimeters <- sim$firePerimeters[sim$firePerimeters$YEAR > 2000 &
                                               sim$firePerimeters$YEAR < 2012,]
  }

  ## Species equivalencies table -------------------------------------------
  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")

    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)
    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## check spp column to use
    if (P(sim)$sppEquivCol == "Boreal") {
      message(paste("There is no 'sppEquiv' table supplied;",
                    "will attempt to use species listed under 'Boreal'",
                    "in the 'LandR::sppEquivalencies_CA' table"))
    } else {
      if (grepl(P(sim)$sppEquivCol, names(sim$sppEquiv))) {
        message(paste("There is no 'sppEquiv' table supplied,",
                      "will attempt to use species listed under", P(sim)$sppEquivCol,
                      "in the 'LandR::sppEquivalencies_CA' table"))
      } else {
        stop("You changed 'sppEquivCol' without providing 'sppEquiv',",
             "and the column name can't be found in the default table ('LandR::sppEquivalencies_CA').",
             "Please provide conforming 'sppEquivCol', 'sppEquiv' and 'sppColorVect'")
      }
    }

    ## remove empty lines/NAs
    sim$sppEquiv <- sim$sppEquiv[!"", on = P(sim)$sppEquivCol]
    sim$sppEquiv <- na.omit(sim$sppEquiv, P(sim)$sppEquivCol)

    ## add default colors for species used in model
    sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                                  newVals = "Mixed", palette = "Accent")
  } else {
    if (is.null(sim$sppColorVect))
      stop("If you provide 'sppEquiv' you MUST also provide 'sppColorVect'")
  }

  ## Species raster layers -------------------------------------------
  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")

    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)
    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## check spp column to use
    if (P(sim)$sppEquivCol == "Boreal") {
      message(paste("There is no 'sppEquiv' table supplied;",
                    "will attempt to use species listed under 'Boreal'",
                    "in the 'LandR::sppEquivalencies_CA' table"))
    } else {
      if (grepl(P(sim)$sppEquivCol, names(sim$sppEquiv))) {
        message(paste("There is no 'sppEquiv' table supplied,",
                      "will attempt to use species listed under", P(sim)$sppEquivCol,
                      "in the 'LandR::sppEquivalencies_CA' table"))
      } else {
        stop("You changed 'sppEquivCol' without providing 'sppEquiv',",
             "and the column name can't be found in the default table ('LandR::sppEquivalencies_CA').",
             "Please provide conforming 'sppEquivCol', 'sppEquiv' and 'sppColorVect'")
      }
    }

    ## remove empty lines/NAs
    sim$sppEquiv <- sim$sppEquiv[!"", on = P(sim)$sppEquivCol]
    sim$sppEquiv <- na.omit(sim$sppEquiv, P(sim)$sppEquivCol)

    ## add default colors for species used in model
    sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                                  newVals = "Mixed", palette = "Accent")
  } else {
    if (is.null(sim$sppColorVect))
      stop("If you provide 'sppEquiv' you MUST also provide 'sppColorVect'")
  }

  if (!suppliedElsewhere("speciesLayersStart", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    sim$speciesLayersStart <- Cache(loadkNNSpeciesLayers,
                                    dPath = dPath,
                                    rasterToMatch = sim$rasterToMatch,
                                    studyArea = sim$studyArea,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMapStart, LCC.. etc
                                    sppEquiv = sim$sppEquiv,
                                    knnNamesCol = "KNN",
                                    sppEquivCol = P(sim)$sppEquivCol,
                                    thresh = P(sim)$coverThresh,
                                    url = extractURL("speciesLayersStart"),
                                    userTags = c(cacheTags, "speciesLayersStart"),
                                    omitArgs = c("userTags"))
  }

  ## filter sppEquiv according to available layers
  sppEquiv <- which(sim$sppEquiv[[P(sim)$sppEquivCol]] %in% names(sim$speciesLayersStart))
  sppEquiv <- sim$sppEquiv[sppEquiv]

  if (!suppliedElsewhere("speciesLayersEnd", sim)) {
    sim$speciesLayersEnd <- Cache(loadkNNSpeciesLayersValidation,
                                  dPath = dPath,
                                  rasterToMatch = sim$rasterToMatch,
                                  studyArea = sim$studyArea,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMapStart, LCC.. etc
                                  sppEquiv = sppEquiv,
                                  knnNamesCol = "KNN",
                                  sppEquivCol = P(sim)$sppEquivCol,
                                  thresh = P(sim)$coverThresh,
                                  url = extractURL("speciesLayersEnd"),
                                  userTags = c(cacheTags, "speciesLayersEnd"),
                                  omitArgs = c("userTags"))
  }

  ## check which species layers common between the two objects
  commonLayers <- intersect(names(sim$speciesLayersEnd),
                            names(sim$speciesLayersStart))
  if (length(commonLayers)) {
    if (any(!names(sim$speciesLayersStart) %in% commonLayers))
      message(red(paste("The following species were used for simulation, but have no validation data:\n",
                        paste(names(sim$speciesLayersStart)[!names(sim$speciesLayersStart) %in% commonLayers],
                              collapse = ", "),
                        "\nThey will be excluded from the validation, but the user should consider whether they should be",
                        "\nin the simulation, or providing other validation layers.")))
  } else {
    stop("There are no common species between 'speciesLayersStart' and 'speciesLayersEnd'.
         Please check these objects and/or they are being produced")
  }

  ## Biomass layers ----------------------------------------------------
  if (!suppliedElsewhere("rawBiomassMapEnd", sim)) {
    rawBiomassValFileName <- "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
    httr::with_config(config = httr::config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
      #necessary for KNN
      sim$rawBiomassMapEnd <- Cache(prepInputs,
                                    targetFile = rawBiomassValFileName,
                                    url = extractURL("rawBiomassMapEnd"),
                                    destinationPath = asPath(dPath),
                                    fun = "raster::raster",
                                    studyArea = sim$studyArea,
                                    rasterToMatch = sim$rasterToMatch,
                                    useSAcrs = FALSE,
                                    method = "bilinear",
                                    datatype = "INT2U",
                                    filename2 = TRUE,
                                    overwrite = TRUE,
                                    userTags = c(cacheTags, "rawBiomassMapEnd"),
                                    omitArgs = c("userTags"))
    })
  }

  if (!suppliedElsewhere("biomassMap", sim)) {
    sim$biomassMap <- sim$rawBiomassMapStart
  }

  ## Age layers ----------------------------------------------------
  if (!suppliedElsewhere("standAgeMapStart", sim)) {
    httr::with_config(config = httr::config(ssl_verifypeer = 0L), {
      sim$standAgeMapStart <- Cache(LandR::prepInputsStandAgeMap,
                                    destinationPath = dPath,
                                    ageURL = extractURL("standAgeMapStart"),
                                    studyArea = raster::aggregate(sim$studyAreaLarge),
                                    rasterToMatch = sim$rasterToMatchLarge,
                                    filename2 = .suffix("standAgeMapStart.tif", paste0("_", P(sim)$.studyAreaName)),
                                    overwrite = TRUE,
                                    fireURL = extractURL("fireURL"),
                                    fireField = "YEAR",
                                    startTime = start(sim),
                                    userTags = c("prepInputsStandAge_rtm", currentModule(sim), cacheTags),
                                    omitArgs = c("destinationPath", "targetFile", "overwrite",
                                                 "alsoExtract", "userTags"))
    })
  }

  if (!suppliedElsewhere("standAgeMapEnd", sim)) {
    standAgeValFileName <- "NFI_MODIS250m_2011_kNN_Structure_Stand_Age_v1.tif"
    httr::with_config(config = httr::config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
      #necessary for KNN
      sim$standAgeMapEnd <- Cache(prepInputs,
                                  targetFile = standAgeValFileName,
                                  url = extractURL("standAgeMapEnd"),
                                  destinationPath = asPath(dPath),
                                  fun = "raster::raster",
                                  studyArea = sim$studyArea,
                                  rasterToMatch = sim$rasterToMatch,
                                  useSAcrs = FALSE,
                                  method = "bilinear",
                                  datatype = "INT2U",
                                  filename2 = TRUE,
                                  overwrite = TRUE,
                                  userTags = c(cacheTags, "standAgeMapEnd"),
                                  omitArgs = c("userTags"))
    })
  }

  ## Cohort data -------------------------------------------
  if (!suppliedElsewhere("allCohortData", sim)) {
    if (!suppliedElsewhere("simulationOutputs", sim)) {
      stop("If not supplying 'allCohort' then you MUST supply 'simulationOutputs'")
    } else {
      cohortDataOutputs <- sim$simulationOutputs[objectName == "cohortData"]

      ## check that the selected years and reps exist in outputs table
      if (!any(is.na(P(sim)$validationReps))) {
        out <- lapply(P(sim)$validationYears, FUN = function(y, cohortDataOutputs, reps) {
          fileNames <- cohortDataOutputs[saveTime == y, file]
          reps <- paste("rep", reps, sep = "")
          reps <- sub("(rep)([[:digit:]])$", "\\10\\2", reps)
          out <- vapply(reps, FUN = function(x) any(grepl(x, fileNames)), FUN.VALUE = logical(1))
          out2 <- which(!out)
          if (length(out2))
            stop(paste("Missing 'cohortData' outputs for year", y, "and rep(s)", paste(out2, collapse = ", ")))
        }, cohortDataOutputs = cohortDataOutputs, reps = P(sim)$validationReps)
      } else {
        out <- lapply(P(sim)$validationYears, FUN = function(y, cohortDataOutputs) {
          fileNames <- cohortDataOutputs[saveTime == y, file]
          if (!length(fileNames))
            stop(paste("Missing 'cohortData' outputs for year", y))
        }, cohortDataOutputs = cohortDataOutputs)
      }

      ##  subset to validation years
      cohortDataOutputs <- cohortDataOutputs[saveTime %in% P(sim)$validationYears]

      ## subset to validation reps and add reps column
      if (!any(is.na(P(sim)$validationReps))) {
        repsStr <- paste("rep", P(sim)$validationReps, sep = "")
        repsStr <- sub("(rep)([[:digit:]])$", "\\10\\2", repsStr)
        repsStr <- paste(repsStr, sep = "", collapse = "|" )
        cohortDataOutputs <- cohortDataOutputs[grepl(repsStr, file)]
        cohortDataOutputs[, rep := sub("(.*rep)([[:digit:]]*)(\\/.*)", "\\2", file)]
      } else {
        cohortDataOutputs[, rep := "1"]
      }

      ## check that files exist
      if (any(!vapply(cohortDataOutputs$file, file.exists, FUN.VALUE = logical(1)))) {
        stop("All/some cohortData files cannot be found.
             This can be because the simulation outputs were saved in different paths from the ones
             passed to Biomass_validationKNN in sim$simulationOutputs.
             Please verify paths in 'simulationOutputs' and correct if necessary")
      }

      sim$allCohortData <- rbindlist(fill = TRUE, use.names = TRUE,
                                     l = apply(cohortDataOutputs, MARGIN = 1, FUN = function(x) {
                                       cohortData <- readRDS(x["file"])
                                       cohortData[, year := as.numeric(x["saveTime"])]
                                       cohortData[, rep := as.numeric(x["rep"])]
                                       return(cohortData)
                                     }))
    }
  }

  if (!suppliedElsewhere("pixelGroupMapStk", sim)) {
    if (!suppliedElsewhere("simulationOutputs", sim)) {
      stop("If not supplying 'allCohort' then you MUST supply 'simulationOutputs'")
    } else {
      pixelGroupMapOutputs <- sim$simulationOutputs[objectName == "pixelGroupMap"]

      ## check that the selected years and reps exist in outputs table
      if (!any(is.na(P(sim)$validationReps))) {
        out <- lapply(P(sim)$validationYears, FUN = function(y, pixelGroupMapOutputs, reps) {
          fileNames <- pixelGroupMapOutputs[saveTime == y, file]
          reps <- paste("rep", reps, sep = "")
          reps <- sub("(rep)([[:digit:]])$", "\\10\\2", reps)
          out <- vapply(reps, FUN = function(x) any(grepl(x, fileNames)), FUN.VALUE = logical(1))
          out2 <- which(!out)
          if (length(out2))
            stop(paste("Missing 'pixelGroupMap' outputs for year", y, "and rep(s)", paste(out2, collapse = ", ")))
        }, pixelGroupMapOutputs = pixelGroupMapOutputs, reps = P(sim)$validationReps)
      } else {
        out <- lapply(P(sim)$validationYears, FUN = function(y, pixelGroupMapOutputs) {
          fileNames <- pixelGroupMapOutputs[saveTime == y, file]
          if (!length(fileNames))
            stop(paste("Missing 'pixelGroupMap' outputs for year", y))
        }, pixelGroupMapOutputs = pixelGroupMapOutputs)
      }

      ##  subset to validation years
      pixelGroupMapOutputs <- pixelGroupMapOutputs[saveTime %in% P(sim)$validationYears]

      ## subset to validation reps and add reps column
      if (!any(is.na(P(sim)$validationReps))) {
        repsStr <- paste("rep", P(sim)$validationReps, sep = "")
        repsStr <- sub("(rep)([[:digit:]])$", "\\10\\2", repsStr)
        repsStr <- paste(repsStr, sep = "", collapse = "|" )
        pixelGroupMapOutputs <- pixelGroupMapOutputs[grepl(repsStr, file)]
        pixelGroupMapOutputs[, rep := sub("(.*rep)([[:digit:]]*)(\\/.*)", "\\2", file)]
      } else {
        pixelGroupMapOutputs[, rep := "1"]
      }


      ## check that files exist
      if (any(!vapply(pixelGroupMapOutputs$file, file.exists, FUN.VALUE = logical(1)))) {
        stop("All/some pixelGroupMap files cannot be found.
             This can be because the simulation outputs were saved in different paths from the ones
             passed to Biomass_validationKNN in sim$simulationOutputs.
             Please verify paths in 'simulationOutputs' and correct if necessary")
      }

      sim$pixelGroupMapStk <- apply(pixelGroupMapOutputs, MARGIN = 1, FUN = function(x) {
        pixelGroupMap <- readRDS(x["file"])
        names(pixelGroupMap) <- paste0("year", as.numeric(x["saveTime"]), "_rep", as.numeric(x["rep"]))
        pixelGroupMap
      })
      sim$pixelGroupMapStk <- stack(sim$pixelGroupMapStk)
    }
  }

  return(invisible(sim))
}
