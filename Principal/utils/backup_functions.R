Load10X_Spatial_2 <- function(
    data.dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    bin.size = NULL,
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = NULL,
    ...
) {
  # if more than one directory is passed in
  if (length(x = data.dir) > 1) {
    # party on with the first value
    data.dir <- data.dir[1]
    # but also raise a warning
    warning(
      paste0(
        "`data.dir` expects a single value but recieved multiple - ",
        "continuing using the first: '",
        data.dir,
        "'."
      ),
      immediate. = TRUE,
    )
  }
  # if the specified directory does not exist
  if (!file.exists(data.dir)) {
    # raise an error
    stop(paste0("No such file or directory: ", "'", data.dir, "'"))
  }
  
  # if `bin.size` is not set but `data.dir` points to a folder with binned data
  if (is.null(bin.size) & file.exists(paste0(data.dir, "/binned_outputs"))) {
    # point `bin.size` to the "standard" set - i.e. everything in the default
    # output except the 8 um binning because it's a memory hog
    bin.size <- c(16, 8)
  }
  # if `bin.size` is specified
  if(!is.null(bin.size)) {
    # convert `bin.size` to a character vector and pad values to three digits
    bin.size.pretty <- paste0(sprintf("%03d", bin.size), "um")
    # point `data.dirs` to the specified binnings
    data.dirs <- paste0(
      data.dir,
      "/binned_outputs/",
      "square_",
      bin.size.pretty
    )
    # suffix assay/slice names with each bin size
    assay.names <- paste0(assay, ".", bin.size.pretty)
    slice.names <- paste0(slice, ".", bin.size.pretty)
  } else {
    # otherwise just hold onto the top-level directory
    data.dirs <- data.dir
    # and keep the assay/slice names unchanged
    assay.names <- assay
    slice.names <- slice
  }
  
  # read in counts matrices from specified h5 files
  counts.paths <- lapply(data.dirs, file.path, filename)
  counts.list <- lapply(counts.paths, Read10X_h5, ...)
  # maybe convert Cell identifiers to uppercase
  if (to.upper) {
    rownames(counts) <- lapply(rownames(counts), toupper)
  }
  
  if (is.null(image)) {
    # read in the corresponding images and coordinate mappings
    image.list <- mapply(
      Read10X_Image_2,
      file.path(data.dirs, "spatial"),
      assay = assay.names,
      slice = slice.names,
      MoreArgs = list(filter.matrix = filter.matrix)
    )
  } else {
    # make sure any passed images are in a vector
    image.list <- c(image)
  }
  
  # check that for each counts matrix there is a corresponding image
  if (length(image.list) != length(counts.list)) {
    stop(
      paste0(
        "The number of images does not match the number of counts matrices. ",
        "Ensure each spatial dataset has a corresponding image."
      )
    )
  }
  
  # for each counts matrix, build a Seurat object
  object.list <- mapply(CreateSeuratObject, counts.list, assay = assay.names)
  # associate each counts matrix with its corresponding image
  object.list <- mapply(
    function(
    .object,
    .image,
    .assay,
    .slice
    ) {
      # align the image's identifiers with the object's
      .image <- .image[Cells(.object)]
      # add the image to the corresponding Seurat instance
      .object[[.slice]] <- .image
      return (.object)
    },
    object.list,
    image.list,
    assay.names,
    slice.names
  )
  # merge the Seurat instances - each assay should have unique Cell identifiers
  object <- merge(
    object.list[[1]],
    y = object.list[-1]
  )
  
  return(object)
}

Read10X_Image_2 <- function(
    image.dir,
    image.name = "tissue_lowres_image.png",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE
) {
  image <- png::readPNG(
    source = file.path(
      image.dir,
      image.name
    )
  )
  
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(image.dir, "scalefactors_json.json")
  )
  
  coordinates <- Read10X_Coordinates(
    filename = Sys.glob(file.path(image.dir, "*tissue_positions*")),
    filter.matrix
  )
  
  coordinates$imagerow <- as.numeric(coordinates$imagerow)
  coordinates$imagecol <- as.numeric(coordinates$imagecol)
  
  fov <- CreateFOV(
    coordinates[, c("imagerow", "imagecol")],
    type = "centroids",
    radius = scale.factors[["spot"]],
    assay = assay,
    key = Key(slice, quiet = TRUE)
  )
  
  visium.fov <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = image,
    scale.factors = scale.factors
  )
  
  return(visium.fov)
}
