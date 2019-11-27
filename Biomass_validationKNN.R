
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
  version = list(SpaDES.core = "0.2.6", reproducible = "0.2.11",
                 LandR = "0.0.2.9007",
                 Biomass_validationKNN = "0.0.1",
                 Biomass_core = "1.3.2", Biomass_speciesData = "1.0.0",
                 Biomass_borealDataPrep = "1.4.0.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_validationKNN.Rmd"),
  reqdPkgs = list("crayon", "reproducible", "raster",
                  "sf"),
  parameters = rbind(
    defineParameter("coverThresh", "integer", "10", NA, NA,
                    paste("The minimum % cover a species needs to have (per pixel) in the study",
                          "area to be considered present. Should be the same as the one used to obtain",
                          "the species cover layers for simualtion set up.")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA,
                    desc = paste("Vector of length = 1, describing the simulation time at",
                                 "which the first plot event should occur. Set to NA to turn plotting off.")),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default")
    ),
  inputObjects = bind_rows(
    expectsInput("firePerimeters", "sf",
                 desc = paste("A map of fire perimeters in the study area that can be used to exclude pixels",
                              "that have been burnt during the validation period. Defaults to the Canadian",
                              "Wildland Fire Information System 1986-2018 National Burned Area Composite,",
                              "subset to fires between 2001 and 2011 (inclusively)."),
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nbac/nbac_1986_to_2018_20191017.zip"),
    expectsInput("rawBiomassMap", "RasterLayer",
                 desc = paste("total biomass raster layer in study area. Defaults to the Canadian Forestry",
                              "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                              "from 2001. See http://tree.pfc.forestry.ca/NFI_MAP_V0_metadata.xls for metadata"),
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("rawBiomassMapValidation", "RasterLayer",
                 desc = paste("total biomass raster layer in study area used for validation. Defaults to the",
                                      "Canadian Forestry Service, National Forest Inventory, kNN-derived total aboveground",
                                      "biomass map from 2011. See",
                                      "https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990"),
                 sourceURL = "http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/"),
    expectsInput("rstLCChange", "RasterLayer",
                 desc = paste("A map of land cover change types in the study area that can be used to exclude pixels",
                              "that have been disturbed during the validation period. Defaults to Canada's forest",
                              "change national map between 1985-2011 (CFS), subset to years 2001-2011 (inclusively).",
                              "See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information."),
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/C2C_change_type_1985_2011.zip"),
    expectsInput("speciesLayers", "RasterStack",
                 desc = paste("cover percentage raster layers by species in Canada species map.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived species cover maps from 2001, using a cover threshold of 10% -",
                              "see http://tree.pfc.forestry.ca/NFI_MAP_V0_metadata.xls for metadata"),
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput("speciesLayersValidation", "RasterStack",
                 desc = paste("cover percentage raster layers by species in Canada species map used for validation.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived species cover maps from 2011 -",
                              "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = "http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/")
  ),
  outputObjects = bind_rows(
    createsOutput("disturbedIDs", "RasterLayer",
                  desc = paste("pixels that have been disturbed by fire or suffered land-cover changes",
                               "during the validation period. These pixels are excluded form the validation.")),
    createsOutput("rawBiomassMapValidation", "RasterLayer",
                  desc = paste("total biomass raster layer in study area used for validation.",
                               "Filtered to exclude pixels that were disturbed during the validation period")),
    createsOutput("speciesLayersValidation", "RasterStack",
                 desc = paste("cover percentage raster layers by species in Canada species map used for validation.",
                              "Filtered to exclude pixels that were disturbed during the validation period"))
  )
))

doEvent.Biomass_validationKNN = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
    },

    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  ## EXCLUDE DISTURBED PIXELS FROM VALIDATION  -----------------------------------
  ## make a template raster with IDs
  rasterToMatchIDs <- sim$rasterToMatch
  rasterToMatchIDs <- setValues(rasterToMatchIDs, values = 1:ncell(rasterToMatchIDs))

  ## get pixels inside fire perimeters
  inFireIDs <- getValues(mask(rasterToMatchIDs, sim$firePerimeters))
  inFireIDs <- inFireIDs[!is.na(inFireIDs)] ## faster than na.omit

  inLCChangeIDs <- rasterToMatchIDs[!is.na(sim$rstLCChange)]

  ## make vector of pixels that are both in fire perimeters and LCChange
  ## export to sim
  sim$disturbedIDs <- union(inFireIDs, inLCChangeIDs)

  ## exclude these pixels from validation layers
  sim$speciesLayersValidation[sim$disturbedIDs] <- NA
  sim$rawBiomassMapValidation[sim$disturbedIDs] <- NA

  ## return some statistics about excluded pixels
  excludedPixStats <- data.table(noPixels = length(sim$disturbedIDs),
                                 landscapePrc = round(length(sim$disturbedIDs)/
                                                        sum(!is.na(getValues(sim$biomassMap))),
                                                      2) * 100)
  message(blue("Pixels disturbed during the validation period will be excluded from validation, representing a loss of:\n",
               excludedPixStats$noPixels, "pixels or", excludedPixStats$landscapePrc, "% of the landscape."))
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

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
      message("There is no rasterToMatch supplied; will attempt to use rawBiomassMap")
    } else {
      stop("rasterToMatch is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (!suppliedElsewhere("rawBiomassMap", sim) || needRTM) {
    sim$rawBiomassMap <- Cache(prepInputs,
                               targetFile = asPath(basename(rawBiomassMapFilename)),
                               archive = asPath(c("kNN-StructureBiomass.tar",
                                                  "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                               url = extractURL("rawBiomassMap"),
                               destinationPath = dPath,
                               studyArea = sim$studyArea,
                               rasterToMatch = if (!needRTM) sim$rasterToMatch else NULL,
                               # maskWithRTM = TRUE,    ## if RTM not supplied no masking happens (is this intended?)
                               maskWithRTM = if (!needRTM) TRUE else FALSE,
                               ## TODO: if RTM is not needed use SA CRS? -> this is not correct
                               # useSAcrs = if (!needRTM) TRUE else FALSE,
                               useSAcrs = FALSE,     ## never use SA CRS
                               method = "bilinear",
                               datatype = "INT2U",
                               filename2 = TRUE, overwrite = TRUE, userTags = cacheTags,
                               omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
  }
  if (needRTM) {
    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMap
    sim$rasterToMatch <- sim$rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatch)
    sim$rasterToMatch[!is.na(RTMvals)] <- 1

    sim$rasterToMatch <- Cache(writeOutputs, sim$rasterToMatch,
                               filename2 = file.path(cachePath(sim), "rasters", "rasterToMatch.tif"),
                               datatype = "INT2U", overwrite = TRUE,
                               userTags = cacheTags,
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
                                 userTags = cacheTags,
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

    rstLCChangeYr <- Cache(prepInputs,
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
    pixKeep <- getValues(sim$rstLCChange) > 0 &
      getValues(rstLCChangeYr) > 100 &
      getValues(rstLCChangeYr) < 112

    sim$rstLCChange[!pixKeep] <- NA
  }

  if (!compareRaster(sim$rstLCChange,
                     sim$rasterToMatch, stopiffalse = FALSE)) {
    stop("'rstLCChange' and 'rasterToMatch' differ in
         their properties. Please check")
  }


  ## Fire perimeter data ---------------------------------------------------

  if (!suppliedElsewhere("firePerimeters", sim)) {
    firePerimetersFile <- "nbac_1986_to_2018_20191017.shp"
    sim$firePerimeters <- Cache(prepInputs,
                                targetFile = firePerimetersFile,
                                alsoExtract = "similar",
                                archive = asPath("nbac_1986_to_2018_20191017.zip"),
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

  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    sim$speciesLayers <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyArea,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                               sppEquiv = sim$sppEquiv,
                               knnNamesCol = "KNN",
                               sppEquivCol = P(sim)$sppEquivCol,
                               thresh = P(sim)$coverThresh,
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"),
                               omitArgs = c("userTags"))
  }

  ## filter sppEquiv according to available layers
  sppEquiv <- which(sim$sppEquiv[[P(sim)$sppEquivCol]] %in% names(sim$speciesLayers))
  sppEquiv <- sim$sppEquiv[sppEquiv]

  if (!suppliedElsewhere("speciesLayersValidation", sim)) {
    sim$speciesLayersValidation <- Cache(loadkNNSpeciesLayersValidation,
                                         dPath = dPath,
                                         rasterToMatch = sim$rasterToMatch,
                                         studyArea = sim$studyArea,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                                         sppEquiv = sppEquiv,
                                         knnNamesCol = "KNN",
                                         sppEquivCol = P(sim)$sppEquivCol,
                                         thresh = P(sim)$coverThresh,
                                         url = extractURL("speciesLayersValidation"),
                                         userTags = c(cacheTags, "speciesLayersValidation"),
                                         omitArgs = c("userTags"))
  }

  ## check which species layers common between the two objects
  commonLayers <- intersect(names(sim$speciesLayersValidation),
                            names(sim$speciesLayers))
  if (length(commonLayers)) {
    if (any(!names(sim$speciesLayers) %in% commonLayers))
      message(red(paste("The following species were used for simulation, but have no validation data:\n",
                        paste(names(sim$speciesLayers)[!names(sim$speciesLayers) %in% commonLayers],
                              collapse = ", "),
                        "\nThey will be excluded from the validation, but the user should consider whether they should be",
                        "\nin the simulation, or providing other validation layers.")))
  } else {
    stop("There are no common species between 'speciesLayers' and 'speciesLayersValidation'.
         Please check these objects and/or they are being produced")
  }

  if (!compareRaster(sim$speciesLayersValidation,
                     sim$rasterToMatch, stopiffalse = FALSE)) {
    stop("'speciesLayersValidation' and 'rasterToMatch' differ in
         their properties. Please check")
  }

  ## Biomass layer ----------------------------------------------------
  if (!suppliedElsewhere("rawBiomassMapValidation", sim)) {
    ## get all online file names
    fileURLs <- getURL(extractURL("rawBiomassMapValidation"), dirlistonly = TRUE)
    fileNames <- getHTMLLinks(fileURLs)
    rawBiomassValFileName <- grep("Biomass_TotalLiveAboveGround.*.tif$", fileNames, value = TRUE)
    rawBiomassValURL <- paste0(extractURL("rawBiomassMapValidation"), rawBiomassValFileName)

    sim$rawBiomassMapValidation <- Cache(prepInputs,
                                         targetFile = asPath(rawBiomassValFileName),
                                         url = rawBiomassValURL,
                                         destinationPath = asPath(dPath),
                                         fun = "raster::raster",
                                         studyArea = sim$studyArea,
                                         rasterToMatch = sim$rasterToMatch,
                                         useSAcrs = FALSE,
                                         method = "bilinear",
                                         datatype = "INT2U",
                                         filename2 = TRUE,
                                         overwrite = TRUE,
                                         userTags = c(cacheTags, "rawBiomassMapValidation"),
                                         omitArgs = c("userTags"))
  }

  ## check rasters
  if (!compareRaster(sim$rawBiomassMapValidation,
                     sim$rasterToMatch, stopiffalse = FALSE)) {
    stop("'rawBiomassMapValidation' and 'rasterToMatch' differ in
         their properties. Please check")
  }

  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
