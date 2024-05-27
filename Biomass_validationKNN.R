# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_validationKNN",
  description = "Validation module for LandR Biomass predictions of forest succession. Based on Canadian Forest Service KNN maps", #"insert module description here",
  keywords = c("validation", "ecological simulation model",
               "forest dynamics", "forest succession", "data", "prediction"),
  authors = c(person("Ceres", "Barros", email = "ceres.barros@ubc.ca", role = c("aut", "cre")),
              person(c("Eliot"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut"))),
  childModules = character(0),
  version = list(Biomass_validationKNN = "0.0.3"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_validationKNN.Rmd"),
  reqdPkgs = list("achubaty/amc", "crayon", "ggplot2", "ggpubr",
                  "mclust", "terra", "RCurl", "scales", "sf", "XML",
                  # "curl", "httr", ## called directly by this module, but pulled in by LandR (Sep 6th 2022).
                  ## Excluded because loading is not necessary (just installation)
                  "PredictiveEcology/LandR@development (>= 1.1.0.9064)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@development (>= 2.0.2)",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9004)",
                  "PredictiveEcology/SpaDES.tools@development"),
  parameters = rbind(
    defineParameter("coverThresh", "integer", "10", NA, NA,
                    desc = paste("The minimum % cover a species needs to have (per pixel) in the study",
                                 "area to be considered present. Should be the same as the one used to obtain",
                                 "the species cover layers for simulation set up.")),
    defineParameter("deciduousCoverDiscount", "numeric", 0.8418911, NA, NA,
                    desc = paste("This was estimated with data from NWT on March 18, 2020 and may or may not be universal.",
                                 "Should be the same as the one used when preparing `cohortData` in the simulation set up.")),
    defineParameter("LCChangeYr", "integer", NULL, 1900, NA,
                    desc = paste("OPTIONAL. An integer or vector of integers of the validation period years, defining which",
                                 "years of land-cover changes (i.e. disturbances) should be excluded.",
                                 "`NULL` by default, which presumes no subsetting based on years is done internally (either",
                                 "the user supplies a pre-filtered `rstLCChange`, or no filtering is desired). If not `NULL`",
                                 "`rstLCChangeYr` is used to filter disturbed pixels within the specified years.",
                                 "See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information.")),
    defineParameter("minCoverThreshold", "numeric", 5, 0, 100,
                    desc = paste("Cover that is equal to or below this number will be omitted from the dataset",
                                 "Should be the same as the one used when preparing `cohortData` in the simulation set up.")),
    defineParameter("obsDeltaAgeB", "logical", TRUE, NA, NA,
                    desc = paste("When TRUE, the observed changes in biomass and age (deltaB, deltaAge) between",
                                 "the two validation years will be plotted as maps and scatterplots")),
    defineParameter("pixelGroupBiomassClass", "numeric", 100, NA, NA,
                    desc = paste("When assigning `pixelGroup` membership, this defines the resolution of biomass that will be",
                                 "considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same",
                                 "Should be the same as the one used when preparing `cohortData` in the simulation set up.")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    "The column in `sim$sppEquiv` data.table to use as a naming convention"),
    defineParameter("validationReps", "integer", 1:10, NA, NA,
                    desc = paste("The simulation repetitions for the validation. Defaults to 1:10. Set to NA if not using repetitions",
                                 "(i.e. only one run)")),
    defineParameter("validationYears", "integer", c(2001, 2011), NA, NA,
                    desc = "The simulation years for the validation. Defaults to 2001 and 2011. Must select two years"),
    defineParameter(".plotInitialTime", "integer", 1L,
                    desc = paste("If NA plotting is off completely (this includes saving).")),
    defineParameter(".plots", "character", default = c("object", "png"),
                    desc = paste("Passed to `types` in Plots (see ?Plots). There are a few plots that are made within this module, if set.",
                                 "Note that plots (or their data) are saved in file.path(outputPath(sim), 'figures').",
                                 "If `NA`, plotting is off completely (this includes plot saving).")),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between save events"),
    defineParameter(".sslVerify", "integer", as.integer(unname(curl::curl_options("^ssl_verifypeer$"))), NA_integer_, NA_integer_,
                    paste("Passed to `httr::config(ssl_verifypeer = P(sim)$.sslVerify)` when downloading KNN",
                          "(NFI) datasets. Set to 0L if necessary to bypass checking the SSL certificate (this",
                          "may be necessary when NFI's website SSL certificate is not correctly configured).")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If `NA`, a hash of `studyArea` will be used."),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default")
  ),
  inputObjects = bind_rows(
    expectsInput("allCohortData", "data.table",
                 desc = paste("All `cohortData` tables saved during the simulation, particularly for the validation years.",
                              "If not supplied, the module will attempt to retrieve them using the 'simulationOutputs' table")),
    expectsInput("biomassMap", "SpatRaster",
                 desc = paste("total biomass raster layer in study area (in g/m^2), filtered for pixels covered by `cohortData`.",
                              "Only used to calculate total no. of pixels being simulated",
                              "If not supplied, will default to `rawBiomassMapStart`")),
    expectsInput("firePerimeters", "sf",
                 desc = paste("A map of fire perimeters in the study area that can be used to exclude pixels",
                              "that have been burnt during the validation period. If burnt pixels are not to be excluded",
                              "Provide an empty `sf` object with the same properties as the default. Defaults to the latest Canadian",
                              "Wildland Fire Information System National Burned Area Composite,",
                              "subset to fires occuring up to last validation year (inclusively). Source URL determined by `fireURL`"),
                 sourceURL = NA),
    expectsInput("fireURL", "character",
                 desc = paste("A URL to a fire database, such as the Canadian National Fire Database,",
                              "that is a zipped shapefile with fire polygons, an attribute (i.e., a column) named 'Year'.",
                              "If supplied (omitted with NULL or NA), this will be used to 'update' age pixels on `standAgeMap`",
                              "with 'time since fire' as derived from this fire polygons map"),
                 sourceURL = "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"),
    expectsInput("pixelGroupMapStk", "RasterStack",
                 desc = paste("A stack of `pixelGroupMap`s saved during the simulation, particularly for the validation years.",
                              "If not supplied, the module will attempt to make it using the 'simulationOutputs' table")),
    expectsInput("rawBiomassMapStart", "SpatRaster",
                 desc = paste("observed total biomass raster layer in study area at the first year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived",
                              "total aboveground biomass map from 2001 (in ton/ha).",
                              "See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata."),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")),
    expectsInput("rawBiomassMapEnd", "SpatRaster",
                 desc = paste("observed total biomass raster layer in study area at the last year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived",
                              "total aboveground biomass map from 2011 (in ton/ha)",
                              "See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
                                    "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")),
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = paste("A raster of the `studyArea` in the same resolution and projection as `rawBiomassMapStart`.",
                              "This is the scale used for all *outputs* for use in the simulation.")),
    expectsInput("rstLCChange", "SpatRaster",
                 desc = paste("A mask-type map of land cover changes in the study area that can be used to exclude pixels",
                              "that have been disturbed during the validation period. If disturbed pixels are not to be excluded",
                              "Provide an empty sf object with the same properties as the default. Defaults to Canada's forest",
                              "change map between 1985-2011 (CFS), filtered for years 2001-2011 (inclusively)",
                              "and all disturbances collapsed (map only has values of 1 and NA). See `P(sim)$LCChangeYr` parameter",
                              "to change the period of disturbances, and",
                              "https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information."),
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/C2C_change_type_1985_2011.zip"),
    expectsInput("rstLCChangeYr", "SpatRaster",
                 desc = paste("An OPTIONAL map of land cover change years in the study area used to exclude pixels that have",
                              "been disturbed during the validation period. It defaults to Canada's forest",
                              "change year national map between 1985-2011 (CFS). If `P(sim)$LCChangeYr` is not `NULL`,",
                              "this layer is used to filted disturbed pixels that fall within the years specified by ",
                              "`P(sim)$LCChangeYr`. If `P(sim)$LCChangeYr` is `NULL` this layer is not used.",
                              "See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information."),
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/C2C_change_year_1985_2011.zip"),
    expectsInput("simulationOutputs", "data.table",
                 desc = paste("An OPTIONAL table listing simulation outputs (as passed to `spades()`, or `experiment`)",
                              "that will be used to make `allCohortData`, `pixelGroupMapStk`,",
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
                 desc = paste("observed percent cover raster layers by species in Canada used for validation",
                              "at the last year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived species cover maps from 2011 -",
                              "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/")),
    expectsInput("sppColorVect", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in sim$sppEquiv[[sim$sppEquivCol]],",
                              "and should also contain a color for 'Mixed'"),
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See `LandR::sppEquivalencies_CA`."),
    expectsInput("standAgeMapStart", "SpatRaster",
                 desc =  paste("observed stand age map in study area, at the first year of the validation period",
                               "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                               "kNN-derived biomass map from 2001 -",
                               "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("standAgeMapEnd", "SpatRaster",
                 desc = paste("observed stand age raster layer in study area, at the last year of the validation period.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory, kNN-derived stand age",
                              "map from 2011. See",
                              "https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
                                    "NFI_MODIS250m_2011_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("studyArea", "SpatVector",
                 desc = paste("Polygon to use as the study area. Must be provided by the user")),
  ),
  outputObjects = bind_rows(
    createsOutput("logLikelihood", "data.table",
                  desc = paste("A table of negative sum log-likelihood values calculated for different variables and averaged",
                               "across repetitions. At the moment, log-likelihood values are calculated for biomass (landscape- and",
                               "pixel-level), species presences and dominance (lanscape-level) and deltaB (landscape- and pixel-level.",
                               "For biomass and count data (presences/dominance, we assume an underlying multinomial distribution,",
                               "and for deltaB a multivariate Gaussian distribution - note that the later is still under development.")),
    createsOutput("landscapeMAD", "data.table",
                  desc = paste("Mean absolute deviance values calculated on landscape-level relative abundances, species presences and dominance,",
                               "and deltaB, per repetition and year (except for deltaB, which is integrated across years)")),
    createsOutput("landscapeVars", "data.table",
                  desc = paste("A table containing observed and simulated landscape-averaged variables used for validation",
                               "(by year and repetition, 'rep', in the case of simulated data), namely: species relative abundances",
                               "('relAbund'), species presenses ('count'), species dominance (as in no. pixels where a given species,",
                               "has higher 'relAbund'; 'countDom') and species changes in biomass, as 2011 minus 2001 ('deltaB').",
                               "Observed data rows are labelled as 'observed' in 'dataType' column. In species dominance, pixels",
                               "with >= 2 species with max(B) and pixels with no B  are classified as 'Mixed' and 'No veg.',",
                               "respectively in the 'speciesCode' column - note that this is 'vegType' column in `pixelCohortData`.")),
    createsOutput("pixelCohortData", "data.table",
                  desc = paste("A table containing observed and simulated pixel-level data (by year and repetition, 'rep',",
                               "in the case of simulated data) on species biomass (summed across cohorts, 'B'),",
                               "total pixel biomass ('pixelB'), average biomass-weighted pixel age ('pixelAge'),",
                               "species relative abundance (calculated as B/pixelB, 'relativeAbund'), species dominance",
                               "(the species with max(B), 'vegType'), and lanscape-wide biomass ('landscapeB').",
                               "Observed data columns are suffixed with 'Obsrvd'. In species dominance, pixels with >= 2" ,
                               "species with max(B) (i.e. 'noDoms' >= 2) are classified as 'Mixed'." )),
    createsOutput("pixelMAD", "data.table",
                  desc = paste("Mean absolute deviance values calculated on pixel-level relative abundances and deltaB,",
                               "per repetition and year (except for deltaB, which is integrated across years)")),
    createsOutput("pixelVars", "data.table",
                  desc = paste("The same as `landscapeVars`, but variables are calculated at the pixel-level")),
    createsOutput("rstDisturbedPix", "SpatRaster",
                  desc = paste("Raster of pixel IDs (as a mask) that have been disturbed by fire or suffered land-cover",
                               "changes during the validation period. These pixels are excluded form the validation.")),
    createsOutput("rawBiomassMapStart", "SpatRaster",
                  desc = paste("observed total biomass raster layer in study area at the first year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("rawBiomassMapEnd", "SpatRaster",
                  desc = paste("observed total biomass raster layer in study area at the last year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("speciesLayersStart", "SpatRaster",
                  desc = paste("observed percent cover raster layers by species in Canada",
                               "at the first year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("speciesLayersEnd", "SpatRaster",
                  desc = paste("observed percent cover raster layers by species in Canada",
                               "at the last year of the validation period.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("standAgeMapStart", "SpatRaster",
                  desc =  paste("observed stand age map in study area, at the first year of the validation period",
                                "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("standAgeMapEnd", "SpatRaster",
                  desc = paste("observed stand age map in study area, at the last year of the validation period",
                               "Filtered to exclude pixels that were disturbed during the validation period"))

  )
))

doEvent.Biomass_validationKNN = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)

      if (any(grepl("screen", P(sim)$.plots))) {
        dev(height = 7, width = 12)   ## don't overwrite other plots, open new window
        mod$statsWindow <- dev.cur()
        if (P(sim)$obsDeltaAgeB) {
          mod$mapWindow <- mod$statsWindow + 1 ## this window should to be made first if need be
          mod$landscapeWindow <- mod$mapWindow + 1
        } else {
          mod$landscapeWindow <- mod$statsWindow + 1
        }

        mod$pixelWindow <- mod$landscapeWindow + 1
      }

      sim <- scheduleEvent(sim, eventType = "calculateValidVars", eventTime = times(sim)$start,
                           eventPriority = 1, moduleName = currentModule(sim))
      sim <- scheduleEvent(sim, eventType = "validationStats", eventTime = times(sim)$start,
                           eventPriority = 2, moduleName = currentModule(sim))

      ## all plotting events have the same priority - they don't depend on each other.
      if (P(sim)$obsDeltaAgeB) {
        sim <- scheduleEvent(sim, eventType = "obsDeltaMaps", eventTime = times(sim)$start,
                             eventPriority = 3, moduleName = currentModule(sim))
      }

      if (anyPlotting(P(sim)$.plots)) {
        sim <- scheduleEvent(sim, eventType = "landscapeWidePlots", eventTime = times(sim)$start,
                             eventPriority = 3, moduleName = currentModule(sim))
        sim <- scheduleEvent(sim, eventType = "pixelLevelPlots", eventTime = times(sim)$start,
                             eventPriority = 3, moduleName = currentModule(sim))
        sim <- scheduleEvent(sim, eventType = "deltaBComparisons", eventTime = times(sim)$start,
                             eventPriority = 3, moduleName = currentModule(sim))
      }
    },
    calculateValidVars = {
      sim <- calculateValidVarsEvent(sim)
    },
    validationStats = {
      sim <- validationStatsEvent(sim)
    },
    obsDeltaMaps = {
      sim <- obsrvdDeltaMapsEvent(sim)
    },
    landscapeWidePlots = {
      sim <- landscapeWidePlotsEvent(sim)
    },
    pixelLevelPlots = {
      sim <- pixelLevelPlotsEvent(sim)
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
  mod$plotPath <- file.path(outputPath(sim), "figures")
  if (!dir.exists(mod$plotPath))
    dir.create(mod$plotPath, recursive = TRUE)

  ## make internal module reps variable
  mod$validationReps <- if (!any(is.na(P(sim)$validationReps)))
    P(sim)$validationReps else
      1L

  ## CHECK RASTER LAYERS AGAINST RASTERTOMATCH -----------------------------------
  if (!.compareRas(sim$biomassMap, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$biomassMap <- postProcess(sim$biomassMap, rasterToMatch = sim$rasterToMatch)
  }

  if (!.compareRas(sim$rawBiomassMapStart, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$rawBiomassMapStart <- postProcess(sim$rawBiomassMapStart, rasterToMatch = sim$rasterToMatch)
  }
  if (!.compareRas(sim$rawBiomassMapEnd, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$rawBiomassMapEnd <- postProcess(sim$rawBiomassMapEnd, rasterToMatch = sim$rasterToMatch)
  }

  if (!.compareRas(sim$speciesLayersStart, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$speciesLayersStart <- postProcess(sim$speciesLayersStart, rasterToMatch = sim$rasterToMatch)
  }

  if (!.compareRas(sim$speciesLayersEnd, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$speciesLayersEnd <- postProcess(sim$speciesLayersEnd, rasterToMatch = sim$rasterToMatch)
  }

  if (!.compareRas(sim$standAgeMapStart, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$standAgeMapStart <- postProcess(sim$standAgeMapStart, rasterToMatch = sim$rasterToMatch)
  }

  if (!.compareRas(sim$standAgeMapEnd, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$standAgeMapEnd <- postProcess(sim$standAgeMapEnd, rasterToMatch = sim$rasterToMatch)
  }

  if (!.compareRas(sim$rstLCChange, sim$rasterToMatch, stopOnError = FALSE)) {
    sim$rstLCChange <- postProcess(sim$rstLCChange, rasterToMatch = sim$rasterToMatch)
  }

  if (!is.null(P(sim)$LCChangeYr)) {
    if (!.compareRas(sim$rstLCChangeYr, sim$rasterToMatch, stopOnError = FALSE)) {
      sim$rstLCChangeYr <- postProcess(sim$rstLCChangeYr, rasterToMatch = sim$rasterToMatch)
    }
  }

  ## EXCLUDE DISTURBED PIXELS FROM VALIDATION  -----------------------------------
  ## make a template raster with IDs
  rasterToMatchIDs <- sim$rasterToMatch
  rasterToMatchIDs[] <- 1:ncell(rasterToMatchIDs)

  ## get pixels inside fire perimeters
  if (!all(st_is_empty(sim$firePerimeters))) {
    inFireIDs <- as.vector(mask(rasterToMatchIDs, sim$firePerimeters)[])
    inFireIDs <- inFireIDs[!is.na(inFireIDs)] ## faster than na.omit
  } else {
    inFireIDs <- integer(0)
  }

  ## only keep pixels that have been disturbed during the validation period
  ## convert years to the map's format
  if (!is.null(P(sim)$LCChangeYr)) {
    yrs <- P(sim)$LCChangeYr - 1900
    pixKeep <- !is.na(as.vector(sim$rstLCChange[])) &
      as.vector(sim$rstLCChangeYr[]) %in% yrs

    sim$rstLCChange[!pixKeep] <- NA
    sim$rstLCChangeYr[!pixKeep] <- NA
  }

  ## get pixels inside LCC pixels
  inLCChangeIDs <- rasterToMatchIDs[!is.na(sim$rstLCChange)]

  ## make vector of pixels that are both in fire perimeters and LCChange
  ## export to sim
  mod$disturbedIDs <- union(inFireIDs, inLCChangeIDs)

  sim$rstDisturbedPix <- sim$rasterToMatch
  sim$rstDisturbedPix[] <- NA
  sim$rstDisturbedPix[mod$disturbedIDs] <- 1

  ## exclude these pixels from validation layers
  sim$speciesLayersEnd[mod$disturbedIDs] <- NA
  sim$rawBiomassMapEnd[mod$disturbedIDs] <- NA
  sim$standAgeMapEnd[mod$disturbedIDs] <- NA

  sim$speciesLayersStart[mod$disturbedIDs] <- NA
  sim$rawBiomassMapStart[mod$disturbedIDs] <- NA
  sim$standAgeMapStart[mod$disturbedIDs] <- NA

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

  validationDataStart <- LandR:::.createCohortData(pixelTable,
                                                   sppColumns = grep("cover\\.", names(pixelTable), value = TRUE),
                                                   rescale = TRUE,
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
                                                 sppColumns = grep("cover\\.", names(pixelTable), value = TRUE),
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

  ## summarize allPixelCohortData to pixel totalB per species and
  ## and biomass-averaged pixel age
  pixelCohortData <- allPixelCohortData[, .(B, sum(B), age),
                                        by = .(rep, year, pixelIndex, speciesCode)]
  pixelCohortData <- pixelCohortData[, pixelAge := sum(age * B, na.rm = TRUE) / sum(B, na.rm = TRUE),
                                     by = .(rep, year, pixelIndex)]
  ## drop unnecessary columns and remove separate cohorts
  pixelCohortData[, B := V2] ## overwrite
  pixelCohortData[, `:=`(V2 = NULL, age = NULL)]
  pixelCohortData <- unique(pixelCohortData)

  ## MERGE OBSERVED AND SIMULATED DATA TABLES  -----------------------------------
  ## add observed data to pixelCohortData
  ## note that some pixelIndex X spp combinations are lacking because the observed data
  ## has spp in some pixels that are not found in the simulation data,
  ## and vice-versa. To make sure that the observed pixel X spp combinations are added to each
  ## rep/year the tables need to be extended - otherwise sometimes the observed data
  ## is only joined to some reps/years, making observed averages "vary" across reps/years
  combinationsStart <- as.data.table(expand.grid(list(speciesCode = unique(pixelCohortData$speciesCode),
                                                      pixelIndex = unique(c(validationDataStart$pixelIndex, pixelCohortData$pixelIndex)),
                                                      rep = mod$validationReps,
                                                      year = P(sim)$validationYears[1])))

  combinationsEnd <- as.data.table(expand.grid(list(speciesCode = unique(pixelCohortData$speciesCode),
                                                    pixelIndex = unique(c(validationDataEnd$pixelIndex, pixelCohortData$pixelIndex)),
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
  sim$validationData <- sim$validationData[pixelIndex %in% pixelCohortData$pixelIndex]

  ## convert NAs to 0s
  cols <- c("cover", "age", "B", "totalBiomass")
  sim$validationData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]

  ## change names before joining.
  setnames(sim$validationData,
           old = c("cover", "age", "B", "totalBiomass"),
           new = c("coverObsrvd", "pixelAgeObsrvd", "BObsrvd", "pixelBObsrvd"))
  ## reorder columns
  cols <- c("rep", "year", "pixelIndex", "speciesCode", "coverObsrvd",
            "pixelAgeObsrvd", "BObsrvd", "pixelBObsrvd")
  sim$validationData <- sim$validationData[, ..cols]

  ## merge and keep all combos
  pixelCohortData <- pixelCohortData[sim$validationData,
                                     on = c("rep", "year", "pixelIndex", "speciesCode")]

  ## remove disturbed pixels
  pixelCohortData <- pixelCohortData[!pixelIndex %in% mod$disturbedIDs]

  ## convert NAs to 0s
  cols <- c("pixelAge", "B")
  pixelCohortData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]

  ## calculate simulated pixelAge, pixelB and relative B
  pixelCohortData[, pixelB := asInteger(sum(B)),
                  .(rep, year, pixelIndex)]
  pixelCohortData[, relativeAbund := B/pixelB]
  pixelCohortData[pixelB == 0, relativeAbund := 0]

  ## replace inserted pixelB/AgeObsrvd 0s (coming from merge) with actual pixel biomass/age
  pixelCohortData[, `:=`(pixelBObsrvd = max(pixelBObsrvd),
                         pixelAgeObsrvd = max(pixelAgeObsrvd)),
                  by = .(rep, year, pixelIndex)]

  ## calculate observed relative B
  pixelCohortData[, relativeAbundObsrvd := BObsrvd/pixelBObsrvd]
  pixelCohortData[pixelBObsrvd == 0, relativeAbundObsrvd := 0]

  ## classify pixels by dominant species (i,e, veg type) - will need to be corrected for compeeting dominants
  ## these match with model outputs, I checked
  ## note2: mixed pixels get a "mixed" type
  pixelCohortData[pixelB > 0, vegType := speciesCode[which.max(relativeAbund)],
                  by = .(year, rep, pixelIndex)]
  pixelCohortData[pixelBObsrvd > 0, vegTypeObsrvd := speciesCode[which.max(relativeAbundObsrvd)],
                  by = .(year, rep, pixelIndex)]


  ## get number of dominant species -- note that pixels with 0 B will appear as mixed initially, but are fixed below
  tempDT <- pixelCohortData[, list(noDoms = sum(relativeAbund == max(relativeAbund)),
                                   noDomsObsrvd = sum(relativeAbundObsrvd == max(relativeAbundObsrvd))),
                            by = .(year, rep, pixelIndex)]

  pixelCohortData <- tempDT[pixelCohortData, on = .(year, rep, pixelIndex)]
  pixelCohortData[pixelB == 0, noDoms := 0]
  pixelCohortData[pixelBObsrvd == 0, noDomsObsrvd := 0]

  pixelCohortData[noDoms > 1, vegType := "Mixed"]
  pixelCohortData[noDomsObsrvd > 1, vegTypeObsrvd := "Mixed"]

  pixelCohortData[is.na(vegType), vegType := "No veg."]
  pixelCohortData[is.na(vegTypeObsrvd), vegTypeObsrvd := "No veg."]

  ## calculate some landscape metrics
  pixelCohortData[, `:=`(landscapeB = sum(B),
                         landscapeBObsrvd = sum(BObsrvd)),
                  by = .(rep, year)]

  ## reorder column names
  cols <- c(grep("Obsrvd", names(pixelCohortData), value = TRUE, invert = TRUE),
            grep("Obsrvd", names(pixelCohortData), value = TRUE))
  pixelCohortData <- pixelCohortData[, ..cols]

  ## assert and export to sim -- vegType cols can have NAs
  assertPixelCohortDataValid(pixelCohortData)
  sim$pixelCohortData <- pixelCohortData

  ## make labels for plots
  mod$speciesLabels <- equivalentName(unique(sim$pixelCohortData$speciesCode), sim$sppEquiv,
                                      column = "EN_generic_short")
  names(mod$speciesLabels) <- unique(sim$pixelCohortData$speciesCode)

  ## EXCLUDE PIXELS WHERE OBSERVED PIXEL B OR PIXEL AGE DECREASED -----------------
  ## Only keep pixels where pixel age AND pixel B increased, or remained the same,
  ## as we cannot account for disturbances that may not have been captured from sat data
  ## or measurement errors that yielded too high B in 2001. It is unlikely that in 10 years we have a large proportion
  ## of pixels seeing reduced pixel age and B due to death from long-age.
  ## THIS HAS BEEN DESACTIVATED FOR NOW.

  if (FALSE) {
    year1 <- P(sim)$validationYears[1]
    year2 <- P(sim)$validationYears[2]
    tempDT <- unique(sim$pixelCohortData[, .(year, rep, pixelIndex, pixelAgeObsrvd, pixelBObsrvd)])
    tempDT <- tempDT[, .(pixelDeltaAgeObsrvd = pixelAgeObsrvd[year == year2] - pixelAgeObsrvd[year == year1],
                         pixelDeltaBObsrvd = pixelBObsrvd[year == year2] - pixelBObsrvd[year == year1]),
                     by = .(rep, pixelIndex)]
    pixToKeep <- unique(tempDT[pixelDeltaBObsrvd >= 0 & pixelDeltaAgeObsrvd >= 0,
                               pixelIndex])
  } else {
    pixToKeep <- unique(sim$pixelCohortData$pixelIndex)
  }

  ## return some statistics about excluded pixels
  pixToRm <- unique(c(setdiff(unique(sim$pixelCohortData$pixelIndex), pixToKeep),
                      mod$disturbedIDs))
  excludedPixStats <- data.table(noPixels = length(pixToRm),
                                 landscapePrc = round(length(pixToRm) /
                                                        sum(!is.na(as.vector(sim$biomassMap[]))),
                                                      2) * 100)
  message(blue("Pixels disturbed during the validation period will be excluded from validation,\n",
               "representing a loss of:", excludedPixStats$noPixels, "pixels or",
               excludedPixStats$landscapePrc, "% of the initial simulated landscape."))

  ## keep aforementioned pixels only
  sim$pixelCohortData <- sim$pixelCohortData[pixelIndex %in% pixToKeep]

  ## clean up and free memory
  rm(pixelTable, pixelCohortData, combinationsStart, combinationsEnd, validationDataStart, validationDataEnd)
  .gc()

  return(invisible(sim))
}

calculateValidVarsEvent <- function(sim)  {
  assertPixelCohortDataValid(sim$pixelCohortData)

  ## LANDSCAPE-LEVEL COMPARISONS --------------------
  ## calculate landscape-wide summary variables:
  ## relative species abundances
  ## species presences (count)
  ## dominant species presences (countDom)
  ## delta biomass
  ## relative abundances and no. pixels with a species
  landscapeVars <- sim$pixelCohortData
  landscapeVars <- landscapeVars[, list(B = sum(B),
                                        BObsrvd = sum(BObsrvd),
                                        landscapeB = unique(landscapeB),
                                        landscapeBObsrvd = unique(landscapeBObsrvd),
                                        landRelativeAbund = sum(B)/unique(landscapeB),
                                        landRelativeAbundObsrvd = sum(BObsrvd)/unique(landscapeBObsrvd),
                                        count = sum(B > 0),
                                        countObsrvd = sum(BObsrvd > 0)),
                                 by = .(rep, year, speciesCode)]
  landscapeVars <- melt.data.table(landscapeVars, measure.vars = list(c("B", "BObsrvd"),
                                                                      c("landRelativeAbund", "landRelativeAbundObsrvd"),
                                                                      c("count", "countObsrvd"),
                                                                      c("landscapeB", "landscapeBObsrvd")),
                                   variable.name = "dataType",
                                   value.name = c("B", "relAbund", "count", "landscapeB"))
  landscapeVars[, dataType := ifelse(dataType == "1", "simulated", "observed")]

  ## no. pixels with a certain dominant species
  ## note: don't use melt, because dominant spp differ between valid and simul data.
  ## simulated and observed differ in no. of pixels in year 1...
  ## this is because B is adjusted using a statistical model
  tempData <- unique(sim$pixelCohortData[, .(year, rep, pixelIndex, vegType)])
  tempData[, dataType := "simulated"]
  tempData2 <- unique(sim$pixelCohortData[, .(year, rep, pixelIndex, vegTypeObsrvd)])
  tempData2[, dataType := "observed"]
  setnames(tempData2, "vegTypeObsrvd", "vegType")

  tempData3 <- rbind(tempData, tempData2)
  tempData3 <- tempData3[, list(countDom = length(pixelIndex)),
                         by = .(rep, year, dataType, vegType)]

  landscapeVars <- landscapeVars[tempData3, on = c("rep", "year", "dataType", "speciesCode==vegType")]
  rm(tempData, tempData2, tempData3)

  ## now make sure that all species/year combinations exist for both datasets
  ## and add zeroes to species that had no biomass/presences
  combos <- expand.grid(dataType = unique(landscapeVars$dataType),
                        rep = unique(landscapeVars$rep),
                        year = unique(landscapeVars$year),
                        speciesCode = unique(landscapeVars$speciesCode)) %>%
    as.data.table(.)
  landscapeVars <- landscapeVars[combos, on = .(dataType, rep, year, speciesCode)]
  cols <- c("B", "relAbund", "count", "countDom")
  landscapeVars <- landscapeVars[, (cols) := replaceNAs(.SD, val = 0), .SDcols = cols]

  landscapeVars[, landscapeB := unique(landscapeB[!is.na(landscapeB)]), by = .(dataType, rep, year)]

  ## delta biomass per species and across the lanscape
  year1 <- P(sim)$validationYears[1]
  year2 <- P(sim)$validationYears[2]

  ## across landscape (calculate landscape-wide changes in total B per species/pixel)
  tempData <- landscapeVars[, list(deltaB = as.numeric(B[which(year == year2)] - B[which(year == year1)]),
                                   landDeltaB = as.numeric(landscapeB[which(year == year2)] - landscapeB[which(year == year1)])),
                            by = .(dataType, rep, speciesCode)]
  ## put "landscape" into speciesCode
  cols <- setdiff(names(tempData), "deltaB")
  tempData2 <- tempData[, ..cols]
  setnames(tempData2, "landDeltaB", "deltaB")
  tempData2 <- unique(tempData2[, speciesCode := "landscape"])

  cols <- setdiff(names(tempData), "landDeltaB")
  tempData <- rbind(tempData[, ..cols], tempData2, use.names = TRUE)

  ## rbind as there is no year correspondence
  landscapeVars <- rbind(landscapeVars, tempData, use.names = TRUE, fill = TRUE)
  rm(tempData, tempData2)


  ## PIXEL-LEVEL COMPARISONS --------------------
  ## relative species abundances
  ## delta biomass
  pixelVars <- sim$pixelCohortData
  pixelVars <- pixelVars[, .(rep, year, pixelIndex, speciesCode, B, BObsrvd,
                             pixelB, pixelBObsrvd,
                             relativeAbund, relativeAbundObsrvd)]
  pixelVars <- melt.data.table(pixelVars,
                               measure.vars = list(c("B", "BObsrvd"),
                                                   c("pixelB", "pixelBObsrvd"),
                                                   c("relativeAbund", "relativeAbundObsrvd")),
                               variable.name = "dataType", value.name = c("B", "pixelB", "relAbund"))
  pixelVars[, dataType := ifelse(dataType == "1", "simulated", "observed")]

  ## now make sure that all species/year/pixel combinations exist for both datasets
  ## and add zeroes to species that had no biomass/presences
  combos <- expand.grid(dataType = unique(pixelVars$dataType),
                        rep = unique(pixelVars$rep),
                        year = unique(pixelVars$year),
                        pixelIndex = unique(pixelVars$pixelIndex),
                        speciesCode = unique(pixelVars$speciesCode)) %>%
    as.data.table(.)
  pixelVars <- pixelVars[combos, on = .(dataType, rep, year, pixelIndex, speciesCode)]
  cols <- c("B", "relAbund", "pixelB")
  pixelVars <- pixelVars[, (cols) := replaceNAs(.SD, val = 0), .SDcols = cols]

  pixelVars[, pixelB := unique(pixelB[!is.na(pixelB)]), by = .(dataType, rep, year, pixelIndex)]

  ## delta biomass per species and pixel (i.e. stand)
  tempData <- pixelVars[, list(deltaB = as.numeric(B[which(year == year2)] - B[which(year == year1)]),
                               pixelDeltaB = as.numeric(pixelB[which(year == year2)] - pixelB[which(year == year1)])),
                        by = .(dataType, rep, pixelIndex, speciesCode)]

  ## put "pixel" into speciesCode
  cols <- setdiff(names(tempData), "deltaB")
  tempData2 <- tempData[, ..cols]
  setnames(tempData2, "pixelDeltaB", "deltaB")
  tempData2 <- unique(tempData2[, speciesCode := "pixel"])

  cols <- setdiff(names(tempData), "pixelDeltaB")
  tempData <- rbind(tempData[, ..cols], tempData2, use.names = TRUE)

  ## rbind as there is no year correspondence
  pixelVars <- rbind(pixelVars, tempData, use.names = TRUE, fill = TRUE)
  rm(tempData, tempData2)

  ## export to sim
  sim$landscapeVars <- landscapeVars
  sim$pixelVars <- pixelVars

  return(invisible(sim))
}

validationStatsEvent <- function(sim) {
  ## MEAN ABSOLUTE DEVIATION -------------------
  ## for each rep, calculate the mean absolute deviation of each pixel to the observed average
  ## note that "observed average" is calculated across all pixels

  ## PIXEL-LEVEL
  ## calculate observed means first
  pixelMAD <- sim$pixelVars[dataType == "observed",
                            list(meanRelAbundObs = mean(relAbund),
                                 meanDeltaBObs = mean(deltaB)),
                            by = .(year, speciesCode)]

  ## checks
  if (NROW(pixelMAD[is.na(year) & !is.na(meanRelAbundObs)])) {
    stop("There should be no relAbund values for 'NA' years (i.e. the deltaB rows)")
  }
  if (NROW(pixelMAD[!is.na(year) & is.na(meanRelAbundObs)])) {
    stop("There are relAbund NAs in non-deltaB rows")
  }

  ## join simulated data and compute MAD
  pixelMAD <- pixelMAD[sim$pixelVars[dataType == "simulated"], on = .(year, speciesCode)]
  pixelMAD <- pixelMAD[, list(meanAbsDevRelAbund = 1/.N * sum(abs(relAbund - meanRelAbundObs)),
                              meanAbsDevDeltaB = 1/.N * sum(abs(deltaB - meanDeltaBObs))),
                       by = .(rep, year, pixelIndex, speciesCode)]

  ## LANDSCAPE-LEVEL
  ## for each rep, calculate the mean absolute deviation of the simulated
  ## landscape-level variable to the observed value
  ## note that "observed average" here is not an actual average
  ## as there is only one observed landscape value
  landscapeMAD <- sim$landscapeVars[dataType == "observed",
                                    list(meanRelAbundObs = unique(relAbund),
                                         meanCountObs = unique(count),
                                         meanCountDomObs = unique(countDom),
                                         meanDeltaBObs = unique(deltaB)),
                                    by = .(year, speciesCode)]
  ## checks
  if (any(duplicated(landscapeMAD[, .(year, speciesCode)]))) {
    stop("There are non-unique observed values for relAbund/count/countDom/deltaB for year X speciesCode combinations")
  }

  ## join simulated data and compute MAD
  landscapeMAD <- landscapeMAD[sim$landscapeVars[dataType == "simulated"], on = .(year, speciesCode)]
  landscapeMAD <- landscapeMAD[, list(meanAbsDevRelAbund = 1/.N * sum(abs(relAbund - meanRelAbundObs)),
                                      meanAbsDevCount = 1/.N * sum(abs(count - meanCountObs)),
                                      meanAbsDevCountDom = 1/.N * sum(abs(countDom - meanCountDomObs)),
                                      meanAbsDevDeltaB = 1/.N * sum(abs(deltaB - meanDeltaBObs))),
                               by = .(rep, year, speciesCode)]


  ## PLOTS
  ## melt and add all labels to factor for equal colours
  plotData <- melt(pixelMAD, measure.vars = c("meanAbsDevRelAbund", "meanAbsDevDeltaB"),
                   value.name = "MAD")
  plotData$variable <- factor(plotData$variable,
                              levels = c("meanAbsDevRelAbund", "meanAbsDevCount", "meanAbsDevCountDom", "meanAbsDevDeltaB"),
                              labels = c("frac('species B', 'total/pixel B')",
                                         "paste('No. of pixels')",
                                         "paste('No. of pixels ')",
                                         "g/m^2"))

  colLabels <- list("frac('species B', 'total/pixel B')" = "rel. abundance",
                    "paste('No. of pixels')" = "presences",
                    "paste('No. of pixels ')" = "dominance",
                    "g/m^2" = bquote(paste(Delta, B)))

  Plots(data = plotData, fn = MADplots,
        filename = "pixelMAD", path = file.path(mod$plotPath),
        deviceArgs = list(width = 8, height = 10, units = "in", res = 300),
        xvar = "speciesCode", yvar = "MAD", colourvar = "variable",
        xlabs = mod$speciesLabels, collabs = colLabels)

  plotData <- melt(landscapeMAD,
                   measure.vars = c("meanAbsDevRelAbund", "meanAbsDevCount", "meanAbsDevCountDom", "meanAbsDevDeltaB"),
                   value.name = "MAD")
  plotData$variable <- factor(plotData$variable,
                              levels = unique(plotData$variable),
                              labels = c("frac('species B', 'total/pixel B')",
                                         "paste('No. of pixels')",
                                         "paste('No. of pixels ')",
                                         "g/m^2"))

  Plots(data = plotData, fn = MADplots,
        filename = "landscapeMAD", path = file.path(mod$plotPath),
        deviceArgs = list(width = 8, height = 10, units = "in", res = 300),
        xvar = "speciesCode", yvar = "MAD", colourvar = "variable",
        xlabs = mod$speciesLabels, collabs = colLabels)


  ## LOG-LIKELIHOOD -------------------
  ## Biomass at landscape-level
  spBObs <- na.omit(sim$landscapeVars[dataType == "observed", .(year, rep, speciesCode, B)])
  spBSim <- na.omit(sim$landscapeVars[dataType == "simulated", .(year, rep, speciesCode, B)])
  spBObs <- dcast(spBObs, ... ~ speciesCode, value.var = "B", fill = 0)
  spBSim <- dcast(spBSim, ... ~ speciesCode, value.var = "B", fill = 0)

  ## probs == 0 lead to infinite values. give them v. small values
  cols <- intersect(names(spBSim), names(mod$speciesLabels))
  spBSim[, (cols) := lapply(.SD, function(x) {x[x == 0] <- 1E-6; x}), .SDcols = cols]

  landscapeLogLikB <- NegSumLogLikWrapper(spBObs, spBSim, reps = P(sim)$validationReps,
                                          years = P(sim)$validationYears, cols = cols,
                                          varType = "biomass")

  ## Biomass at pixel-level
  spBObs <- na.omit(sim$pixelVars[dataType == "observed", .(year, rep, pixelIndex, speciesCode, B)])
  spBSim <- na.omit(sim$pixelVars[dataType == "simulated", .(year, rep, pixelIndex, speciesCode, B)])
  spBObs <- dcast(spBObs, ... ~ speciesCode, value.var = "B", fill = 0)
  spBSim <- dcast(spBSim, ... ~ speciesCode, value.var = "B", fill = 0)

  ## probs == 0 lead to infinite values. give them v. small values
  cols <- intersect(names(spBSim), names(mod$speciesLabels))
  spBSim[, (cols) := lapply(.SD, function(x) {x[x == 0] <- 1E-6; x}), .SDcols = cols]

  pixelLogLikB <- NegSumLogLikWrapper(spBObs, spBSim, reps = P(sim)$validationReps,
                                      years = P(sim)$validationYears, cols = cols,
                                      varType = "biomass")

  ## Counts at landscape-level
  spCountObs <- na.omit(sim$landscapeVars[dataType == "observed", .(year, rep, speciesCode, count)])
  spCountSim <- na.omit(sim$landscapeVars[dataType == "simulated", .(year, rep, speciesCode, count)])
  spCountObs <- dcast(spCountObs, ... ~ speciesCode, value.var = "count", fill = 0)
  spCountSim <- dcast(spCountSim, ... ~ speciesCode, value.var = "count", fill = 0)

  ## probs == 0 lead to infinite values. give them v. small values
  cols <- intersect(names(spCountSim), names(mod$speciesLabels))
  spCountSim[, (cols) := lapply(.SD, function(x) {x[x == 0] <- 1E-6; x}), .SDcols = cols]

  landscapeLogLikCount <- NegSumLogLikWrapper(spCountObs, spCountSim, reps = P(sim)$validationReps,
                                              years = P(sim)$validationYears, cols = cols,
                                              varType = "counts")

  ## Dominance at landscape-level
  spCountObs <- na.omit(sim$landscapeVars[dataType == "observed", .(year, rep, speciesCode, countDom)])
  spCountSim <- na.omit(sim$landscapeVars[dataType == "simulated", .(year, rep, speciesCode, countDom)])
  spCountObs <- dcast(spCountObs, ... ~ speciesCode, value.var = "countDom", fill = 0)
  spCountSim <- dcast(spCountSim, ... ~ speciesCode, value.var = "countDom", fill = 0)

  ## probs == 0 lead to infinite values. give them v. small values
  cols <- intersect(names(spCountSim), c(names(mod$speciesLabels), "Mixed", "No veg."))
  spCountSim[, (cols) := lapply(.SD, function(x) {x[x == 0] <- 1E-6; x}), .SDcols = cols]

  landscapeLogLikCountDom <- NegSumLogLikWrapper(spCountObs, spCountSim, reps = P(sim)$validationReps,
                                                 years = P(sim)$validationYears, cols = cols,
                                                 varType = "counts")


  ## deltaB at landscape-level
  message(red("Computation of log-likelihood for deltaB is still under development"))
  spDeltaObs <- na.omit(sim$landscapeVars[dataType == "observed", .(rep, speciesCode, deltaB)])
  spDeltaSim <- na.omit(sim$landscapeVars[dataType == "simulated", .(rep, speciesCode, deltaB)])
  spDeltaObs <- dcast(spDeltaObs, ... ~ speciesCode, value.var = "deltaB")
  spDeltaSim <- dcast(spDeltaSim, ... ~ speciesCode, value.var = "deltaB")

  cols <- intersect(names(spDeltaSim), c(names(mod$speciesLabels), "pixel"))

  landscapeLogLikDeltaB <- NegSumLogLikWrapper(spDeltaObs, spDeltaSim,
                                               reps = P(sim)$validationReps, cols = cols,
                                               varType = "delta")
  colnames(landscapeLogLikDeltaB) <- paste0(P(sim)$validationYears[2],"-", P(sim)$validationYears[1])

  ## deltaB at pixel-level
  spDeltaObs <- na.omit(sim$pixelVars[dataType == "observed", .(rep, pixelIndex, speciesCode, deltaB)])
  spDeltaSim <- na.omit(sim$pixelVars[dataType == "simulated", .(rep, pixelIndex, speciesCode, deltaB)])
  spDeltaObs <- dcast(spDeltaObs, ... ~ speciesCode, value.var = "deltaB")
  spDeltaSim <- dcast(spDeltaSim, ... ~ speciesCode, value.var = "deltaB")

  cols <- intersect(names(spDeltaSim), c(names(mod$speciesLabels), "pixel"))

  pixelLogLikDeltaB <- NegSumLogLikWrapper(spDeltaObs, spDeltaSim,
                                           reps = P(sim)$validationReps, cols = cols,
                                           varType = "delta")
  colnames(pixelLogLikDeltaB) <- paste0(P(sim)$validationYears[2],"-", P(sim)$validationYears[1])


  ## average across reps and make data.table
  landscapeLikelihood <- rbindlist(
    list("biomass" = data.table(landscapeLogLikB)[, lapply(.SD, mean)],
         "presences" = data.table(landscapeLogLikCount)[, lapply(.SD, mean)],
         "dominance" = data.table(landscapeLogLikCountDom)[, lapply(.SD, mean)],
         "deltaB" = data.table(landscapeLogLikDeltaB)[, lapply(.SD, mean)]),
    use.names = TRUE, idcol = "variable", fill = TRUE)

  pixelLikelihood <- rbindlist(
    list("biomass" = data.table(pixelLogLikB)[, lapply(.SD, mean)],
         "deltaB" = data.table(pixelLogLikDeltaB)[, lapply(.SD, mean)]),
    use.names = TRUE, idcol = "variable", fill = TRUE)

  logLikelihood <- rbindlist(list("landscape" = landscapeLikelihood, "pixel" = pixelLikelihood),
                             use.names = TRUE, idcol = "variable", fill = TRUE)

  ## export to sim
  sim$pixelMAD <- pixelMAD
  sim$landscapeMAD <- landscapeMAD
  sim$logLikelihood <- logLikelihood

  return(invisible(sim))
}

obsrvdDeltaMapsEvent <- function(sim) {
  ## MAPS OF OBSERVED CHANGES IN PIXEL B AND AGE - RAW DATA -------------------
  pixelDeltaBObsrvdRas <- (sim$rawBiomassMapEnd - sim$rawBiomassMapStart) * 100 ## to rescale to to/ha
  pixToNA <- setdiff(1:ncell(pixelDeltaBObsrvdRas), sim$pixelCohortData$pixelIndex)
  pixelDeltaBObsrvdRas[pixToNA] <- NA

  pixelDeltaAgeObsrvdRas <- sim$standAgeMapEnd - sim$standAgeMapStart
  pixelDeltaAgeObsrvdRas[pixToNA] <- NA

  ## what is the relationship between the two?
  pixelDeltaObsrvdData <- na.omit(data.table(pixelIndex = 1:ncell(pixelDeltaBObsrvdRas),
                                             pixelDeltaBObsrvd = as.vector(pixelDeltaBObsrvdRas[]),
                                             pixelDeltaAgeObsrvd = as.vector(pixelDeltaAgeObsrvdRas[])))

  plot1 <- ggplot(pixelDeltaObsrvdData,
                  aes(x = pixelDeltaAgeObsrvd, y = pixelDeltaBObsrvd)) +
    geom_point() +
    # stat_smooth(method = "lm") +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    labs(y = expression(paste("observed ", Delta, "B")))

  plot2 <- ggplot(pixelDeltaObsrvdData,
                  aes(y = pixelDeltaBObsrvd)) +
    geom_boxplot() +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE) +
    theme(axis.title.y = element_blank(), axis.text = element_blank())

  plot3 <- ggplot(pixelDeltaObsrvdData,
                  aes(y = pixelDeltaAgeObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste("observed ", Delta, "age"))) +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE) +
    theme(axis.text.y = element_blank())

  plot4 <- ggarrange(plot1, plot2,  plot3, nrow = 2, ncol = 2, align = "hv",
                     widths = c(1, 0.5), heights = c(1, 0.5))

  ## delta biomass for the "supposed" age increment - all over the place
  year1 <- P(sim)$validationYears[1]
  year2 <- P(sim)$validationYears[2]
  yearGap <- year2 - year1
  plot5 <- ggplot(pixelDeltaObsrvdData[pixelDeltaAgeObsrvd == yearGap],
                  aes(y = pixelDeltaBObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste("observed ", Delta, "B")),
         title = bquote(atop("observed" ~ Delta ~ "B", "in pixels that aged" ~ .(yearGap)))) +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE)

  ## OBSERVED CHANGES IN PIXEL B AND AGE - AFTER ADJUSTMENTS -------
  ## by adjustments we mean, the data cleanup that we replicate from Biomass_borealDataPrep
  plotData <- unique(sim$pixelCohortData[, .(year, pixelIndex, pixelBObsrvd, pixelAgeObsrvd)])
  plotData <- plotData[, list(pixelDeltaBObsrvd = unique(pixelBObsrvd[which(year == year2)]) - unique(pixelBObsrvd[which(year == year1)]),
                              pixelDeltaAgeObsrvd = unique(pixelAgeObsrvd[which(year == year2)]) - unique(pixelAgeObsrvd[which(year == year1)])),
                       , by = pixelIndex]

  if (any(duplicated(plotData$pixelIndex)))
    stop("There should not be duplicated pixels in observed data")

  pixelDeltaBObsrvdAdj[] <- NA
  pixelDeltaBObsrvdAdj[plotData[, pixelIndex]] <- plotData[, pixelDeltaBObsrvd]

  pixelDeltaAgeObsrvdAdj[] <- NA
  pixelDeltaAgeObsrvdAdj[plotData[, pixelIndex]] <- plotData[, pixelDeltaAgeObsrvd]


  plot6 <- ggplot(plotData,
                  aes(x = pixelDeltaAgeObsrvd, y = pixelDeltaBObsrvd)) +
    geom_point() +
    # stat_smooth(method = "lm") +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    labs(y = expression(paste("observed ", Delta, "B", " - adjusted")))

  plot7 <- ggplot(pixelDeltaObsrvdData,
                  aes(y = pixelDeltaBObsrvd)) +
    geom_boxplot() +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE) +
    theme(axis.title.y = element_blank(), axis.text = element_blank())

  plot8 <- ggplot(pixelDeltaObsrvdData,
                  aes(y = pixelDeltaAgeObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste("observed ", Delta, "age", " - adjusted"))) +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE) +
    theme(axis.text.y = element_blank())

  plot9 <- ggarrange(plot6, plot7,  plot8, nrow = 2, ncol = 2, align = "hv",
                     widths = c(1, 0.5), heights = c(1, 0.5))

  ## delta biomass for the "supposed" age increment - all over the place
  plot10 <- ggplot(plotData[pixelDeltaAgeObsrvd == yearGap],
                   aes(y = pixelDeltaBObsrvd)) +
    geom_boxplot() +
    coord_flip() +
    plotTheme(base_size = 12, margin = FALSE, majorYlines = FALSE) +
    labs(y = expression(paste("observed ", Delta, "B", " - adjusted")),
         title = bquote(atop("observed" ~ Delta ~ "B - adjusted", "in pixels that aged" ~ .(yearGap))))

  ## stack rasters for plotting
  plotStack <- stack(pixelDeltaBObsrvdRas, pixelDeltaAgeObsrvdRas,
                     pixelDeltaBObsrvdAdj, pixelDeltaAgeObsrvdAdj)
  names(plotStack) <- c("pixel_deltaB", "pixel_deltaAge",
                        "pixel_deltaB_adjusted", "pixel_deltaAge_adjusted")

  if (anyPlotting(P(sim)$.plots)) {
    if (any(grepl("screen", P(sim)$.plots))) {
      dev(mod$mapWindow)
      clearPlot()
    }

    Plots(plotStack,
          new = TRUE, fn = quickPlot::Plot,
          filename = "deltaB_Age_Maps", path = file.path(mod$plotPath),
          deviceArgs = list(width = 7, height = 7, units = "in", res = 300))

    if (any(grepl("screen", P(sim)$.plots))) {
      dev(mod$statsWindow)
      clearPlot()
      plotDeltaStats <- ggarrange(plot4, plot5, plot9, plot10)
      Plots(plotDeltaStats, types = "screen", new = TRUE, title = "Delta_stats") ## save as separate graphs bellow
    }

    noScreenTypes <- setdiff(P(sim)$.plots, "screen")
    if (length(noScreenTypes)) {
      Plots(plot4, types = noScreenTypes, fn = quickPlot::Plot,
            filename = "observedDeltaBDeltaAge",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 7, height = 5, units = "in", res = 300))
      Plots(plot5, fn = quickPlot::Plot,
            types = noScreenTypes, filename = "observedDeltaB_yearGap",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 5, height = 4, units = "in", res = 300))
      Plots(plot9, fn = quickPlot::Plot,
            types = noScreenTypes, filename = "observedDeltaBDeltaAge_ADJ",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 7, height = 5, units = "in", res = 300))
      Plots(plot10, fn = quickPlot::Plot,
            types = noScreenTypes, filename = "observedDeltaB_yearGapADJ",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 5, height = 4, units = "in", res = 300))
    }
  }

  return(invisible(sim))
}

landscapeWidePlotsEvent <- function(sim) {
  ## Landscape-wide relative abundances
  cols <- c("speciesCode", "relAbund", "dataType", "year")
  plot11 <- ggplot(data = na.omit(sim$landscapeVars[dataType == "simulated", ..cols]),
                   aes(x = speciesCode, y = relAbund)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    stat_summary(data = na.omit(sim$landscapeVars[dataType == "observed", ..cols]),
                 aes(x = speciesCode, y = relAbund, colour = "observed"),
                 fun = "mean", geom = "point", size = 2) +
    scale_x_discrete(labels = mod$speciesLabels) +
    scale_color_manual(values = c("observed" = "red3")) +
    plotTheme(base_size = 12, margin = FALSE, legend = "bottom", x.text.angle = 45) +
    facet_wrap(~ year) +
    labs(title = "Species relative abundances",
         x = "", y = expression(over("species B", "total B")),
         colour = "")

  ## no. pixels with a species
  cols <- c("speciesCode", "count", "dataType", "year")
  plot12 <- ggplot(data = na.omit(sim$landscapeVars[dataType == "simulated", ..cols]),
                   aes(x = speciesCode, y = count)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    stat_summary(data = na.omit(sim$landscapeVars[dataType == "observed", ..cols]),
                 aes(x = speciesCode, y = count, colour = "observed"),
                 fun = "mean", geom = "point", size = 2) +
    scale_x_discrete(labels = mod$speciesLabels) +
    scale_color_manual(values = c("observed" = "red3")) +
    plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
    facet_wrap(~ year) +
    labs(title = "Species presences", x = "", y = "No. of pixels",
         colour = "", fill = "")

  ## no. pixels with a certain dominant species
  cols <- c("speciesCode", "countDom", "dataType", "year")
  plot13 <- ggplot(data = na.omit(sim$landscapeVars[dataType == "simulated", ..cols]),
                   aes(x = speciesCode, y = countDom)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    geom_point(data = na.omit(sim$landscapeVars[dataType == "observed", ..cols]),
               aes(x = speciesCode, y = countDom, colour = "observed"), size = 2) +
    scale_x_discrete(labels = mod$speciesLabels) +
    scale_color_manual(values = c("observed" = "red3")) +
    plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
    facet_wrap(~ year) +
    labs(title = "Dominant species' presences",
         x = "", y = "No. of pixels", fill = "", colour = "")

  maxPixels <- sum(!is.na(as.vector(sim$biomassMap[])))
  plotLandscapeComp <- ggarrange(plot11 + scale_y_continuous(limits = c(0,1)),
                                 plot12 + scale_y_continuous(limits = c(0, maxPixels)),
                                 plot13 + scale_y_continuous(limits = c(0, maxPixels)),
                                 common.legend = TRUE, legend = "bottom",
                                 nrow = 2, ncol = 2)

  if (anyPlotting(P(sim)$.plots)) {
    if (any(grepl("screen", P(sim)$.plots))) {
      dev(mod$landscapeWindow)
      clearPlot()
      Plots(plotLandscapeComp, types = "screen",
            fn = quickPlot::Plot,
            title = "Landscape_averaged_comparisons", new = TRUE)
    }

    noScreenTypes <- setdiff(P(sim)$.plots, "screen")
    if (length(noScreenTypes)) {
      Plots(plot11, fn = quickPlot::Plot,
            types = noScreenTypes, filename = "LandscapeComparisons_relB",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 7, height = 5, units = "in", res = 300))
      Plots(plot12, fn = quickPlot::Plot,
            types = noScreenTypes, filename = "LandscapeComparisons_PresAbs",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 7, height = 5, units = "in", res = 300))
      Plots(plot13, fn = quickPlot::Plot,
            types = noScreenTypes, filename = "LandscapeComparisons_Dom",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 7, height = 5, units = "in", res = 300))
    }
  }

  return(invisible(sim))
}

pixelLevelPlotsEvent <- function(sim) {
  ## Pixel-level relative abundances
  cols <- c("speciesCode", "relAbund", "dataType", "year")
  plot14 <- ggplot(data = na.omit(sim$pixelVars[, ..cols]),
                   aes(x = speciesCode, y = relAbund, fill = dataType )) +
    geom_boxplot() +
    scale_x_discrete(labels = mod$speciesLabels) +
    plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
    facet_wrap(~ year) +
    labs(title = "Species relative abundances", fill = "",
         x = "", y = expression(over("species B", "pixel B")))


  if (anyPlotting(P(sim)$.plots)) {
    if (any(grepl("screen", P(sim)$.plots))) {
      dev(mod$pixelWindow)
      clearPlot()
    }
    Plots(plot14, fn = quickPlot::Plot,
          title = "Pixel_level_comparisons", new = TRUE,
          filename = "PixelComparisons_relB",
          path = file.path(mod$plotPath),
          deviceArgs = list(width = 10, height = 5, units = "in", res = 300))
  }

  return(invisible(sim))
}

deltaBComparisonsEvent <- function(sim) {
  ## COMPARISONS OF DELTA B ------------------
  ## landscape-wide deltaB
  cols <- c("speciesCode", "deltaB", "dataType", "rep")
  plot15 <- ggplot(data = na.omit(sim$landscapeVars[dataType == "simulated", ..cols]),
                   aes(x = speciesCode, y = deltaB, group = rep)) +
    stat_summary(fun = "mean", geom = "bar") +
    stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
    geom_point(data = na.omit(sim$landscapeVars[dataType == "observed", ..cols]),
               aes(x = speciesCode, y = deltaB, colour = "observed"), size = 2) +
    scale_x_discrete(labels = mod$speciesLabels) +
    scale_color_manual(values = c("observed" = "red3")) +
    plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
    labs(title = "Landscape-level", colour = "",
         x = "", y = expression(paste(Delta, "B")))

  ## pixel-level deltaB
  cols <- c("speciesCode", "deltaB", "dataType", "rep")
  plot16 <- ggplot(data = na.omit(sim$pixelVars[, ..cols]),
                   aes(x = speciesCode, y = deltaB, fill = dataType)) +
    geom_boxplot(aes(alpha = speciesCode == "pixel")) +
    scale_x_discrete(labels = c(mod$speciesLabels, "pixel" = "pixel")) +
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 1.0), guide = "none") +
    plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
    labs(title = "Pixel-level", fill = "",
         x = "", y = expression(paste(Delta, "B")))

  simObsDeltaBPlot <- ggarrange(plot15, plot16 + labs(y = " \n "),
                                ncol = 2)

  if (anyPlotting(P(sim)$.plots)) {
    if (any(grepl("screen", P(sim)$.plots))) {
      dev.set(mod$statsWindow)
      Plots(simObsDeltaBPlot, fn = quickPlot::Plot,
            types = "screen", title = "observedDelta", new = TRUE)
    }

    noScreenTypes <- setdiff(P(sim)$.plots, "screen")
    if (length(noScreenTypes)) {
      Plots(plot15, fn = quickPlot::Plot,
            types = noScreenTypes,
            title = "observedDeltaLandscape", new = TRUE,
            filename = "LandscapeComparisons_deltaB",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 5, height = 5, units = "in", res = 300))
      Plots(plot16, fn = quickPlot::Plot,
            types = noScreenTypes,
            title = "observedDeltaPixel", new = TRUE,
            filename = "PixelComparisons_deltaB",
            path = file.path(mod$plotPath),
            deviceArgs = list(width = 5, height = 5, units = "in", res = 300))
    }
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
    stop("Please provide a 'studyArea' polygon")
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- reproducible::studyAreaName(sim$studyArea)
    message("The .studyAreaName is not supplied; derived name from sim$studyArea: ",
            params(sim)[[currentModule(sim)]][[".studyAreaName"]])
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

  ## Biomass layers ----------------------------------------------------
  if (!suppliedElsewhere("rawBiomassMapStart", sim) ||
      (is.null(sim$rawBiomassMapStart) && needRTM)) { ## needs to be in sim now for RTM
    rawBiomassMapFilename <- "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
    # httr::with_config(config = httr::config(ssl_verifypeer =  P(sim)$.sslVerify), {
    #necessary for KNN
    sim$rawBiomassMapStart <- prepRawBiomassMap(
      targetFile = rawBiomassMapFilename,
      url = extractURL("rawBiomassMapStart"),
      studyAreaName = P(sim)$.studyAreaName,
      cacheTags = cacheTags,
      to = if (!needRTM) sim$rasterToMatch else sim$studyArea,
      projectTo = if (!needRTM) NULL else NA, ## don't project to SA if RTMs not present
      destinationPath = dPath,
      filename2 = .suffix("rawBiomassMapStart.tif", paste0("_", P(sim)$.studyAreaName)),
      userTags = c(cacheTags, "rawBiomassMapStart"))
    # })
  }

  if (!suppliedElsewhere("rawBiomassMapEnd", sim)) {
    rawBiomassValFileName <- "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
    # httr::with_config(config = httr::config(ssl_verifypeer =  P(sim)$.sslVerify), {
      #necessary for KNN
      sim$rawBiomassMapEnd <- prepRawBiomassMap(
        targetFile = rawBiomassValFileName,
        url = extractURL("rawBiomassMapEnd"),
        studyAreaName = P(sim)$.studyAreaName,
        cacheTags = cacheTags,
        to = if (!needRTM) sim$rasterToMatch else sim$studyArea,
        projectTo = if (!needRTM) NULL else NA, ## don't project to SA if RTMs not present
        destinationPath = dPath,
        filename2 = .suffix("rawBiomassMapEnd.tif", paste0("_", P(sim)$.studyAreaName)),
        userTags = c(cacheTags, "rawBiomassMapEnd"))
    # })
    }

  if (!suppliedElsewhere("biomassMap", sim)) {
    sim$biomassMap <- sim$rawBiomassMapStart
  }

  ## rasterToMatch --------------------------------
  if (needRTM) {
    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMapStart
    sim$rasterToMatch <- sim$rawBiomassMapStart
    RTMvals <- as.vector(sim$rasterToMatch[])
    sim$rasterToMatch[!is.na(RTMvals)] <- 1L

    sim$rasterToMatch <- Cache(writeTo, sim$rasterToMatch,
                               filename2 = .suffix(file.path(dPath, "rasterToMatch.tif"),
                                                   paste0("_", P(sim)$.studyAreaName)),
                               datatype = "INT2U", overwrite = TRUE,
                               userTags = c(cacheTags, "rasterToMatch"),
                               omitArgs = c("userTags"))

  }

  if (!terra::same.crs(sim$studyArea, sim$rasterToMatch)) {
  # if (!identical(crs(sim$studyArea), crs(sim$rasterToMatch))) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- projectTo(sim$studyArea, crs(sim$rasterToMatch))
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
    sim$rstLCChange <- Cache(prepInputs,
                             targetFile = LCChangeFilename,
                             archive = file.path(dPath, "C2C_change_type_1985_2011.zip"),
                             url = extractURL("rstLCChange"),
                             destinationPath = dPath,
                             to = sim$rasterToMatch,
                             cropTo = sim$studyArea,
                             method = "ngb",
                             datatype = "INT2U",
                             filename2 = .suffix("rstLCChange.tif", paste0("_", P(sim)$.studyAreaName)),
                             overwrite = TRUE,
                             userTags = c("rstLCChange", cacheTags),
                             omitArgs = c("destinationPath", "targetFile", "userTags"))
    ## convert to mask
    sim$rstLCChange[!is.na(sim$rstLCChange[])] <- 1
  }

  ## Check that rstLCChange is a mask and matches RTM
  assertRstLCChange(sim$rstLCChange, sim$rasterToMatch)
  if (!suppliedElsewhere("rstLCChangeYr", sim) & !is.null(P(sim)$LCChangeYr)) {
    LCChangeYrFilename <- "C2C_change_year_1985_2011.tif"
    sim$rstLCChangeYr <- Cache(prepInputs,
                               targetFile = LCChangeYrFilename,
                               archive = asPath("C2C_change_year_1985_2011.zip"),
                               url = extractURL("rstLCChangeYr"),
                               destinationPath = dPath,
                               to = sim$rasterToMatch,
                               cropTo = sim$studyArea,
                               method = "ngb",
                               datatype = "INT2U",
                               filename2 = .suffix("rstLCChangeYr.tif", paste0("_", P(sim)$.studyAreaName)),
                               overwrite = TRUE,
                               userTags = c("rstLCChangeYr", cacheTags),
                               omitArgs = c("destinationPath", "targetFile", "userTags"))
  }

  ## Fire perimeter data ---------------------------------------------------

  if (!suppliedElsewhere("firePerimeters", sim)) {
    sim$firePerimeters <- Cache(prepInputs,
                                url = extractURL("fireURL"),
                                destinationPath = dPath,
                                to = sim$rasterToMatch,
                                cropTo = sim$studyArea,
                                datatype = "INT2U",
                                filename2 = .suffix("firePerimeters.shp", paste0("_", P(sim)$.studyAreaName)),
                                overwrite = TRUE,
                                userTags = c("firePerimeters", cacheTags),
                                omitArgs = c("destinationPath", "targetFile", "userTags"))
    ## convert to sf
    sim$firePerimeters <- st_as_sf(sim$firePerimeters)

    ## exclude fire years outside validation period
    sim$firePerimeters <- sim$firePerimeters[sim$firePerimeters$YEAR >= 1986 &
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
    sim$speciesLayersStart <- Cache(prepSpeciesLayers_KNN,
                                    destinationPath = dPath,
                                    outputPath = dPath,
                                    rasterToMatch = sim$rasterToMatch,
                                    studyArea = sim$studyArea,
                                    studyAreaName = P(sim)$.studyAreaName,
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
    sim$speciesLayersEnd <- Cache(prepSpeciesLayers_KNN,
                                  destinationPath = dPath,
                                  outputPath = dPath,
                                  rasterToMatch = sim$rasterToMatch,
                                  studyArea = sim$studyArea,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMapStart, LCC.. etc
                                  studyAreaName = P(sim)$.studyAreaName,
                                  sppEquiv = sppEquiv,
                                  knnNamesCol = "KNN",
                                  sppEquivCol = P(sim)$sppEquivCol,
                                  thresh = P(sim)$coverThresh,
                                  url = extractURL("speciesLayersEnd"),
                                  year = 2011,
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

  ## Age layers ----------------------------------------------------
  if (!suppliedElsewhere("standAgeMapStart", sim)) {
    # httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      sim$standAgeMapStart <- Cache(LandR::prepInputsStandAgeMap,
                                    ageFun = "terra::rast",
                                    destinationPath = dPath,
                                    ageURL = extractURL("standAgeMapStart"),
                                    to = sim$rasterToMatch,
                                    cropTo = aggregate(sim$studyArea),
                                    filename2 = .suffix("standAgeMapStart.tif", paste0("_", P(sim)$.studyAreaName)),
                                    overwrite = TRUE,
                                    fireURL = extractURL("fireURL"),
                                    fireField = "YEAR",
                                    startTime = start(sim),
                                    userTags = c("prepInputsStandAge_rtm", currentModule(sim), cacheTags),
                                    omitArgs = c("destinationPath", "targetFile", "overwrite",
                                                 "alsoExtract", "userTags"))
    # })
  }

  if (!suppliedElsewhere("standAgeMapEnd", sim)) {
    standAgeValFileName <- "NFI_MODIS250m_2011_kNN_Structure_Stand_Age_v1.tif"
    # httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      #necessary for KNN
      sim$standAgeMapEnd <- Cache(prepInputs,
                                  targetFile = standAgeValFileName,
                                  url = extractURL("standAgeMapEnd"),
                                  destinationPath = asPath(dPath),
                                  fun = "raster::raster",
                                  to = sim$rasterToMatch,
                                  cropTo = sim$studyArea,
                                  method = "bilinear",
                                  datatype = "INT2U",
                                  filename2 = .suffix("standAgeMapEnd.tif", paste0("_", P(sim)$.studyAreaName)),
                                  overwrite = TRUE,
                                  userTags = c(cacheTags, "standAgeMapEnd"),
                                  omitArgs = c("userTags"))
    # })
  }

  ## Cohort data -------------------------------------------
  if (!suppliedElsewhere("allCohortData", sim)) {
    if (!suppliedElsewhere("simulationOutputs", sim)) {
      stop("If not supplying 'allCohortData' then you MUST supply 'simulationOutputs'")
    } else {
      cohortDataOutputs <- sim$simulationOutputs[objectName == "cohortData"]

      ## check that the selected years and reps exist in outputs table
      ## make a reps vector that aligns with file naming
      if (!any(is.na(P(sim)$validationReps))) {
        if (all(grepl("rep[[:digit:]]", cohortDataOutputs$file))) {
          reps <- paste("rep", P(sim)$validationReps, sep = "")
        } else {
          stop(paste("Simulation output paths should contain 'rep' in the",
                     "folder or file name to identify different replicate outputs"))
        }

        repsWzeros <- sub("(rep)([[:digit:]])$", "\\10\\2", reps)
        if (all(grepl(paste(repsWzeros, collapse = "|"), cohortDataOutputs$file))) {
          reps <- repsWzeros
        }
      }

      if (!any(is.na(P(sim)$validationReps))) {
        out <- lapply(P(sim)$validationYears, FUN = function(y, cohortDataOutputs, reps) {
          fileNames <- cohortDataOutputs[saveTime == y, file]
          out <- vapply(reps, FUN = function(x) any(grepl(x, fileNames)), FUN.VALUE = logical(1))
          out2 <- which(!out)
          if (length(out2))
            stop(paste("Missing `cohortData` outputs for year", y, "and rep(s)", paste(out2, collapse = ", ")))
        }, cohortDataOutputs = cohortDataOutputs, reps = reps)
      } else {
        out <- lapply(P(sim)$validationYears, FUN = function(y, cohortDataOutputs) {
          fileNames <- cohortDataOutputs[saveTime == y, file]
          if (!length(fileNames))
            stop(paste("Missing `cohortData` outputs for year", y))
        }, cohortDataOutputs = cohortDataOutputs)
      }

      ##  subset to validation years
      cohortDataOutputs <- cohortDataOutputs[saveTime %in% P(sim)$validationYears]

      ## subset to validation reps and add reps column
      if (!any(is.na(P(sim)$validationReps))) {
        repsStr <- paste(reps, sep = "", collapse = "|" )
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
          out <- vapply(reps, FUN = function(x) any(grepl(x, fileNames)), FUN.VALUE = logical(1))
          out2 <- which(!out)
          if (length(out2))
            stop(paste("Missing 'pixelGroupMap' outputs for year", y, "and rep(s)", paste(out2, collapse = ", ")))
        }, pixelGroupMapOutputs = pixelGroupMapOutputs, reps = reps)
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
        repsStr <- paste(reps, sep = "", collapse = "|" )
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
      isSpat <- vapply(sim$pixelGroupMapStk, is, class2 = "SpatRaster", FUN.VALUE = logical(1))
      if (!all(isSpat)) {
        sim$pixelGroupMapStk[!isSpat] <- lapply(sim$pixelGroupMapStk[!isSpat], rast)
      }
      sim$pixelGroupMapStk <- rast(sim$pixelGroupMapStk)  ## stack
    }
  }

  return(invisible(sim))
}
