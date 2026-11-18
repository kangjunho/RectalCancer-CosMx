## ============================================================
## CRC CosMx final single-script workflow
## Figures 1-4 and Supplementary Figures 1-4
## No external R scripts are required.
## ============================================================

options(timeout = 1200, download.file.method = "libcurl")
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = 1200)
options(crc.cosmx.figure4.only = FALSE)

## ============================================================
## PART 1. Figure 1 and Supplementary Figure 1
## ============================================================

## ============================================================
## Figure 1 and Supplementary Figure 1
## Post-radiotherapy locally advanced rectal adenocarcinoma CosMx
##
## Figure 1
## A  Spatial cell-type maps and matched CellComposite images
## B  Cell counts by region
## C  Cell-type proportions by region
## D  Regional log2 fold change in cell-type proportions
##
## Supplementary Figure 1
## A  Marker-expression DotPlot supporting the final annotations
## ============================================================

options(timeout = 1200, download.file.method = "libcurl")
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = 1200)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(20261012)

## Locate assets distributed next to this script. This works when the
## downloaded bundle is sourced from RStudio or with source()/Rscript.
get_current_script_directory <- function() {
  command_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep(
    "^--file=",
    command_arguments,
    value = TRUE
  )

  if (length(file_argument) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_argument[[1L]]),
      mustWork = FALSE
    )))
  }

  frame_files <- vapply(
    sys.frames(),
    function(current_frame) {
      if (!is.null(current_frame$ofile)) {
        as.character(current_frame$ofile)
      } else {
        NA_character_
      }
    },
    character(1)
  )

  frame_files <- frame_files[
    !is.na(frame_files) & nzchar(frame_files)
  ]

  if (length(frame_files) > 0L) {
    return(dirname(normalizePath(
      frame_files[[length(frame_files)]],
      mustWork = FALSE
    )))
  }

  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()
  ) {
    source_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(error_condition) ""
    )

    if (nzchar(source_path)) {
      return(dirname(normalizePath(
        source_path,
        mustWork = FALSE
      )))
    }
  }

  normalizePath(getwd(), mustWork = FALSE)
}

script_directory <- get_current_script_directory()

normalize_sample_key <- function(sample_values) {
  toupper(
    gsub(
      "[^A-Za-z0-9]",
      "",
      trimws(as.character(sample_values))
    )
  )
}

normalize_fov_id <- function(fov_values) {
  fov_text <- trimws(as.character(fov_values))
  normalized_fov <- rep(NA_character_, length(fov_text))

  direct_numeric <- suppressWarnings(as.numeric(fov_text))
  direct_index <- is.finite(direct_numeric)

  normalized_fov[direct_index] <- as.character(
    as.integer(round(direct_numeric[direct_index]))
  )

  text_index <- !direct_index & nzchar(fov_text)

  if (any(text_index)) {
    trailing_digits <- sub(
      "^.*?([0-9]+)[^0-9]*$",
      "\\1",
      fov_text[text_index],
      perl = TRUE
    )

    trailing_numeric <- suppressWarnings(
      as.integer(trailing_digits)
    )

    normalized_fov[text_index] <- ifelse(
      is.na(trailing_numeric),
      NA_character_,
      as.character(trailing_numeric)
    )
  }

  normalized_fov
}

normalize_region_label <- function(region_values) {
  region_text <- tolower(trimws(as.character(region_values)))
  normalized_region <- rep(NA_character_, length(region_text))

  normalized_region[
    grepl("non[ -]?response", region_text)
  ] <- "Non-response region"

  normalized_region[
    is.na(normalized_region) &
      grepl("response", region_text)
  ] <- "Response region"

  normalized_region
}

## ------------------------------------------------------------
## 1. Paths and study settings
## ------------------------------------------------------------

project_dir <- normalizePath(Sys.getenv("CRC_COSMX_PROJECT_DIR", unset = getwd()), mustWork = TRUE)

object_dir <- file.path(project_dir, "6 Analysis", "Data objects")

figure1_dir <- file.path(
  project_dir,
  "6 Analysis",
  "Figures",
  "Figure1"
)

supplementary1_dir <- file.path(
  project_dir,
  "6 Analysis",
  "Figures",
  "Supplementary Figure1"
)

table1_dir <- file.path(
  project_dir,
  "6 Analysis",
  "Tables",
  "Figure1"
)

supplementary_table1_dir <- file.path(
  project_dir,
  "6 Analysis",
  "Tables",
  "Supplementary Table1"
)

for (current_directory in c(
  figure1_dir,
  supplementary1_dir,
  table1_dir,
  supplementary_table1_dir
)) {
  dir.create(
    current_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

analysis_sample <- Sys.getenv("CRC_COSMX_SAMPLE", unset = "PUBLIC_SAMPLE_ID")
target_fovs <- c("4", "6")
pixel_size_um <- 0.12

region_mapping <- c(
  "4" = "Non-response region",
  "6" = "Response region"
)

region_levels <- c(
  "Non-response region",
  "Response region"
)

region_short_labels <- c(
  "Non-response region" = "Non-response",
  "Response region" = "Response"
)

region_colors <- c(
  "Non-response region" = "#B2182B",
  "Response region" = "#2166AC"
)

celltype_levels <- c(
  "Cancer",
  "CAF_apCAF",
  "CAF_ECM",
  "CAF_iCAF",
  "CAF_myCAF",
  "Imm_Tcell",
  "Imm_Cytotoxic_NK/T",
  "Imm_NK",
  "Imm_Myeloid",
  "Imm_Plasma",
  "Endothelial",
  "Pericyte",
  "Unresolved"
)

display_celltype_levels <- setdiff(
  celltype_levels,
  "Unresolved"
)

## Final manuscript palette. Keep this identical in Figures 1-4.
celltype_colors <- c(
  "Cancer" = "#E3261C",
  "CAF_apCAF" = "#239B6E",
  "CAF_ECM" = "#6FC5AC",
  "CAF_iCAF" = "#53AC45",
  "CAF_myCAF" = "#28BF70",
  "Imm_Tcell" = "#C06BF0",
  "Imm_Cytotoxic_NK/T" = "#367EB5",
  "Imm_NK" = "#81BED8",
  "Imm_Myeloid" = "#5AA2CA",
  "Imm_Plasma" = "#ED79B7",
  "Endothelial" = "#686868",
  "Pericyte" = "#979797",
  "Unresolved" = "#BDBDBD"
)

celltype_colors <- celltype_colors[
  !duplicated(names(celltype_colors))
]

## The bundle includes the exact CellComposite images used in the original
## Figure 1 reference PDF. These take precedence over raw-folder discovery.
figure1_asset_directory <- file.path(
  script_directory,
  "Figure1_assets"
)

bundled_cell_composite_files <- c(
  "4" = file.path(
    figure1_asset_directory,
    "Figure1_CellComposite_Non-response.jpg"
  ),
  "6" = file.path(
    figure1_asset_directory,
    "Figure1_CellComposite_Response.jpg"
  )
)

## Optional manual overrides. Raw images named TME3 CellComposite_F004 and
## TME3 CellComposite_F006 are selected first; bundled images are fallbacks.
cell_composite_files <- c(
  "4" = NA_character_,
  "6" = NA_character_
)

## ------------------------------------------------------------
## 2. Reusable helpers
## ------------------------------------------------------------

theme_manuscript <- theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(
      colour = "black",
      face = "bold",
      size = 8.5
    ),
    axis.title = element_text(
      colour = "black",
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 12,
      margin = margin(b = 3)
    ),
    strip.text = element_text(
      face = "bold",
      colour = "black",
      size = 9
    ),
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(size = 8),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.55
    ),
    plot.margin = margin(5, 6, 5, 6)
  )

save_plot_formats <- function(
  plot_object,
  output_directory,
  base_filename,
  width,
  height,
  save_tiff = TRUE
) {
  ggsave(
    filename = file.path(
      output_directory,
      paste0(base_filename, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(
      output_directory,
      paste0(base_filename, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  if (isTRUE(save_tiff)) {
    ggsave(
      filename = file.path(
        output_directory,
        paste0(base_filename, ".tiff")
      ),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
  }
}

save_panel_collection <- function(
  panel_list,
  output_directory,
  prefix,
  width,
  height
) {
  panel_letters <- LETTERS[seq_along(panel_list)]

  for (panel_index in seq_along(panel_list)) {
    save_plot_formats(
      plot_object = panel_list[[panel_index]],
      output_directory = output_directory,
      base_filename = paste0(prefix, panel_letters[panel_index]),
      width = width[panel_index],
      height = height[panel_index],
      save_tiff = FALSE
    )
  }
}

find_cell_composite_file <- function(
  fov_id,
  manual_file = NA_character_
) {
  bundled_fallback <- bundled_cell_composite_files[[fov_id]]

  return_bundled_fallback <- function() {
    if (
      length(bundled_fallback) == 1L &&
        !is.na(bundled_fallback) &&
        file.exists(bundled_fallback)
    ) {
      return(normalizePath(
        bundled_fallback,
        mustWork = TRUE
      ))
    }

    NA_character_
  }

  if (
    length(manual_file) == 1L &&
      !is.na(manual_file) &&
      file.exists(manual_file)
  ) {
    return(normalizePath(manual_file))
  }

  sample_raw_directory <- file.path(
    project_dir,
    "5 Raw data",
    analysis_sample
  )

  search_directories <- unique(c(
    sample_raw_directory,
    file.path(project_dir, "4 Images"),
    file.path(project_dir, "Images"),
    file.path(project_dir, "5 Raw data")
  ))

  search_directories <- search_directories[
    dir.exists(search_directories)
  ]

  if (length(search_directories) == 0L) {
    return(return_bundled_fallback())
  }

  image_files <- unique(unlist(lapply(
    search_directories,
    function(current_directory) {
      list.files(
        current_directory,
        recursive = TRUE,
        full.names = TRUE,
        pattern = "\\.(png|jpg|jpeg|tif|tiff)$",
        ignore.case = TRUE
      )
    }
  )))

  ## Never reuse images generated by a previous analysis run.
  image_files <- image_files[
    !grepl(
      paste0(
        "6 Analysis",
        .Platform$file.sep,
        "Figures"
      ),
      image_files,
      fixed = TRUE
    )
  ]

  if (length(image_files) == 0L) {
    return(return_bundled_fallback())
  }

  file_basename <- basename(image_files)

  ## Required manuscript source images:
  ## TME3 CellComposite_F004 and TME3 CellComposite_F006.
  exact_tme3_pattern <- paste0(
    "^TME3[ _-]*CellComposite[ _-]*F0*",
    fov_id,
    "\\.(png|jpg|jpeg|tif|tiff)$"
  )

  exact_tme3_files <- image_files[
    grepl(
      exact_tme3_pattern,
      file_basename,
      ignore.case = TRUE,
      perl = TRUE
    )
  ]

  if (length(exact_tme3_files) > 0L) {
    exact_tme3_files <- exact_tme3_files[
      order(
        !grepl(
          analysis_sample,
          exact_tme3_files,
          ignore.case = TRUE
        ),
        nchar(exact_tme3_files),
        exact_tme3_files
      )
    ]

    cat(
      "Using exact TME3 CellComposite for region/FOV ",
      fov_id,
      ":\n",
      exact_tme3_files[[1L]],
      "\n",
      sep = ""
    )

    return(normalizePath(
      exact_tme3_files[[1L]],
      mustWork = TRUE
    ))
  }

  image_search_text <- paste(
    file_basename,
    dirname(image_files)
  )

  composite_match <- grepl(
    "cell.?composite|cellcomposite|composite",
    image_search_text,
    ignore.case = TRUE
  )

  fov_match <- grepl(
    paste0(
      "(fov|f)[^0-9]*0*",
      fov_id,
      "([^0-9]|$)"
    ),
    file_basename,
    ignore.case = TRUE
  ) |
    grepl(
      paste0(
        "(^|[_-])0*",
        fov_id,
        "([_.-]|$)"
      ),
      file_basename,
      ignore.case = TRUE
    )

  matched_files <- image_files[
    composite_match & fov_match
  ]

  if (length(matched_files) == 0L) {
    return(return_bundled_fallback())
  }

  matched_basename <- basename(matched_files)

  candidate_score <-
    100L * grepl(
      analysis_sample,
      matched_files,
      ignore.case = TRUE
    ) +
    50L * grepl(
      "cell.?composite|cellcomposite",
      matched_basename,
      ignore.case = TRUE
    ) +
    25L * grepl(
      paste0(
        "cell.?composite[^0-9]*(fov|f)?[^0-9]*0*",
        fov_id,
        "([^0-9]|$)"
      ),
      matched_basename,
      ignore.case = TRUE
    ) -
    100L * grepl(
      "overlay|segmentation|mask|outline|thumbnail|preview",
      matched_basename,
      ignore.case = TRUE
    )

  matched_files <- matched_files[
    order(
      -candidate_score,
      nchar(matched_files),
      matched_files
    )
  ]

  if (length(matched_files) > 1L) {
    cat(
      "CellComposite candidates for region/FOV ",
      fov_id,
      ":\n",
      paste0("- ", matched_files, collapse = "\n"),
      "\nSelected: ",
      matched_files[[1L]],
      "\n",
      sep = ""
    )
  }

  normalizePath(
    matched_files[[1L]],
    mustWork = TRUE
  )
}

make_placeholder_image_plot <- function(plot_title, missing_path) {
  ggplot() +
    annotate(
      geom = "label",
      x = 0.5,
      y = 0.5,
      label = paste0(
        "CellComposite image not found\n",
        "Set cell_composite_files[\"",
        missing_path,
        "\"] manually"
      ),
      size = 4,
      fontface = "bold"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    labs(title = plot_title) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
}

make_cell_composite_plot <- function(fov_id, plot_title) {
  image_file <- find_cell_composite_file(
    fov_id = fov_id,
    manual_file = cell_composite_files[[fov_id]]
  )

  if (is.na(image_file)) {
    return(
      make_placeholder_image_plot(
        plot_title = plot_title,
        missing_path = fov_id
      )
    )
  }

  cat(
    "CellComposite image for FOV",
    fov_id,
    ":\n",
    image_file,
    "\n"
  )

  image_extension <- tolower(tools::file_ext(image_file))

  if (requireNamespace("magick", quietly = TRUE)) {
    image_raster <- as.raster(
      magick::image_read(image_file)
    )
  } else if (
    image_extension == "png" &&
      requireNamespace("png", quietly = TRUE)
  ) {
    image_raster <- as.raster(
      png::readPNG(image_file)
    )
  } else if (
    image_extension %in% c("jpg", "jpeg") &&
      requireNamespace("jpeg", quietly = TRUE)
  ) {
    image_raster <- as.raster(
      jpeg::readJPEG(image_file)
    )
  } else if (
    image_extension %in% c("tif", "tiff") &&
      requireNamespace("tiff", quietly = TRUE)
  ) {
    image_raster <- as.raster(
      tiff::readTIFF(image_file)
    )
  } else {
    stop(
      "A compatible image reader is required for CellComposite files. ",
      "Install magick, or install png/jpeg/tiff for the corresponding format."
    )
  }

  ggplot() +
    annotation_custom(
      grob = grid::rasterGrob(
        image_raster,
        interpolate = TRUE
      ),
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    labs(title = plot_title) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
}

## ------------------------------------------------------------
## 3. Load the final annotated object
## ------------------------------------------------------------

full_object_file <- file.path(
  object_dir,
  "CRC_CosMx_final_annotated_full_object.rds"
)

if (!file.exists(full_object_file)) stop("Missing required annotated RDS: ", full_object_file)
full_object <- readRDS(full_object_file)

DefaultAssay(full_object) <- "Nanostring"

required_metadata <- c(
  "sample_id",
  "fov",
  "Final_celltype",
  "CenterX_local_px",
  "CenterY_local_px"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(full_object[[]])
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

## ------------------------------------------------------------
## 4. Prepare FOV4/FOV6 metadata
## ------------------------------------------------------------

figure1_data <- as.data.table(
  full_object[[]],
  keep.rownames = "cell"
)

figure1_data[, sample_id_raw := as.character(sample_id)]
figure1_data[, fov_raw := as.character(fov)]
figure1_data[, sample_id := sample_id_raw]
figure1_data[, sample_key := normalize_sample_key(sample_id_raw)]
figure1_data[, fov := normalize_fov_id(fov_raw)]
figure1_data[, Final_celltype := as.character(Final_celltype)]

analysis_sample_key <- normalize_sample_key(analysis_sample)

if ("RT_region" %in% colnames(figure1_data)) {
  figure1_data[
    ,
    Region_from_metadata := normalize_region_label(RT_region)
  ]
} else {
  figure1_data[, Region_from_metadata := NA_character_]
}

figure1_data[
  ,
  Region_from_fov := unname(region_mapping[fov])
]

figure1_data[
  ,
  Region := fifelse(
    !is.na(Region_from_metadata),
    Region_from_metadata,
    Region_from_fov
  )
]

figure1_data[
  ,
  Region := factor(
    Region,
    levels = region_levels
  )
]

cat("\nAvailable sample/FOV labels before standardization:\n")
print(
  unique(
    figure1_data[
      sample_key == analysis_sample_key,
      .(
        sample_id_raw,
        fov_raw,
        normalized_fov = fov,
        RT_region = if (
          "RT_region" %in% colnames(figure1_data)
        ) as.character(RT_region) else NA_character_,
        normalized_region = as.character(Region)
      )
    ]
  )
)

figure1_data <- figure1_data[
  sample_key == analysis_sample_key &
    as.character(Region) %in% region_levels &
    !is.na(CenterX_local_px) &
    !is.na(CenterY_local_px)
]

if (nrow(figure1_data) == 0L) {
  stop(
    paste0(
      "No Figure 1 cells remained after sample/FOV standardization. ",
      "Expected sample key ",
      analysis_sample_key,
      " and RT_region values corresponding to Non-response and Response."
    )
  )
}

figure1_data[
  ,
  x_um := suppressWarnings(
    as.numeric(
      as.character(CenterX_local_px)
    )
  ) * pixel_size_um
]

figure1_data[
  ,
  y_um := suppressWarnings(
    as.numeric(
      as.character(CenterY_local_px)
    )
  ) * pixel_size_um
]

figure1_data <- figure1_data[
  is.finite(x_um) &
    is.finite(y_um)
]

region_cell_check <- figure1_data[
  ,
  .(
    cells = .N,
    min_x_um = min(x_um),
    max_x_um = max(x_um),
    min_y_um = min(y_um),
    max_y_um = max(y_um)
  ),
  by = Region
]

if (!all(region_levels %in% as.character(region_cell_check$Region))) {
  stop(
    "No valid spatial coordinates were retained for region: ",
    paste(
      setdiff(
        region_levels,
        as.character(region_cell_check$Region)
      ),
      collapse = ", "
    )
  )
}

cat("\nFigure 1 spatial cells by region:\n")
print(region_cell_check)

figure1_data[
  !Final_celltype %in% celltype_levels,
  Final_celltype := "Unresolved"
]

figure1_data[
  ,
  Final_celltype := factor(
    Final_celltype,
    levels = celltype_levels
  )
]

## ------------------------------------------------------------
## 5. Figure 1 summary tables
## ------------------------------------------------------------

celltype_count_table <- figure1_data[
  ,
  .(cell_count = .N),
  by = .(
    Region,
    Final_celltype
  )
]

celltype_count_table <- merge(
  CJ(
    Region = factor(region_levels, levels = region_levels),
    Final_celltype = factor(
      celltype_levels,
      levels = celltype_levels
    ),
    unique = TRUE
  ),
  celltype_count_table,
  by = c("Region", "Final_celltype"),
  all.x = TRUE,
  sort = FALSE
)

celltype_count_table[is.na(cell_count), cell_count := 0L]

celltype_count_table[
  ,
  total_cells := sum(cell_count),
  by = Region
]

celltype_count_table[
  ,
  proportion := cell_count / total_cells
]

comparison_table <- dcast(
  celltype_count_table,
  Final_celltype ~ Region,
  value.var = c("cell_count", "proportion"),
  fill = 0
)

comparison_table[
  ,
  log2_response_nonresponse := log2(
    (`proportion_Response region` + 1e-6) /
      (`proportion_Non-response region` + 1e-6)
  )
]

comparison_table[
  ,
  enriched_region := fifelse(
    log2_response_nonresponse >= 0,
    "Response region",
    "Non-response region"
  )
]

comparison_table[
  ,
  Final_celltype := factor(
    Final_celltype,
    levels = celltype_levels
  )
]

fwrite(
  celltype_count_table,
  file.path(
    table1_dir,
    "Figure1_celltype_counts_and_proportions.csv"
  )
)

fwrite(
  comparison_table,
  file.path(
    table1_dir,
    "Figure1_regional_celltype_log2_fold_change.csv"
  )
)

## ------------------------------------------------------------
## 6. Figure 1 panel A: spatial maps and CellComposite images
## ------------------------------------------------------------

make_celltype_spatial_plot <- function(fov_id, plot_title) {
  target_region <- unname(region_mapping[[fov_id]])

  current_data <- copy(
    figure1_data[
      as.character(Region) == target_region
    ]
  )

  if (nrow(current_data) == 0L) {
    stop(
      "No spatial cells found for region ",
      target_region,
      ". Available regions: ",
      paste(
        unique(as.character(figure1_data$Region)),
        collapse = ", "
      )
    )
  }

  ## Pass plain character values to ggplot to avoid factor-level loss in
  ## deeply nested patchwork objects.
  current_data[
    ,
    Plot_celltype := factor(
      as.character(Final_celltype),
      levels = celltype_levels
    )
  ]

  setorder(
    current_data,
    Plot_celltype
  )

  ggplot(
    current_data,
    aes(
      x = x_um,
      y = y_um,
      colour = Plot_celltype
    )
  ) +
    geom_point(
      shape = 16,
      size = 1.65,
      alpha = 1
    ) +
    scale_colour_manual(
      values = celltype_colors,
      limits = celltype_levels,
      breaks = celltype_levels,
      drop = FALSE
    ) +
    coord_fixed(
      ratio = 1,
      xlim = range(current_data$x_um, na.rm = TRUE),
      ylim = range(current_data$y_um, na.rm = TRUE),
      expand = FALSE,
      clip = "off"
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = NULL,
      colour = "Cell type"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12
      ),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
}

p1A_fov4 <- make_celltype_spatial_plot(
  fov_id = "4",
  plot_title = "Non-response region"
) +
  labs(tag = "A")

p1A_fov6 <- make_celltype_spatial_plot(
  fov_id = "6",
  plot_title = "Response region"
)

p1A_composite4 <- make_cell_composite_plot(
  fov_id = "4",
  plot_title = "CellComposite"
)

p1A_composite6 <- make_cell_composite_plot(
  fov_id = "6",
  plot_title = "CellComposite"
)

## Keep the spatial ggplots as direct patchwork children. Wrapping this block
## with wrap_elements() can drop point layers in some patchwork/grid versions.
p1A <- (
  p1A_fov4 | p1A_composite4
) / (
  p1A_fov6 | p1A_composite6
) +
  plot_layout(
    widths = c(1, 1.08),
    heights = c(1, 1)
  )

## ------------------------------------------------------------
## 7. Figure 1 panels B-D
## ------------------------------------------------------------

p1B <- ggplot(
  celltype_count_table[
    as.character(Final_celltype) %in%
      display_celltype_levels
  ],
  aes(
    x = Final_celltype,
    y = cell_count,
    fill = Region
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.72,
    colour = "black",
    linewidth = 0.20
  ) +
  scale_fill_manual(
    values = region_colors,
    labels = region_short_labels,
    drop = FALSE
  ) +
  labs(
    tag = "B",
    title = NULL,
    x = NULL,
    y = "Cell count",
    fill = NULL
  ) +
  theme_manuscript +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 7.5
    )
  )

p1C <- ggplot(
  celltype_count_table[
    as.character(Final_celltype) %in%
      display_celltype_levels
  ],
  aes(
    x = Region,
    y = proportion,
    fill = Final_celltype
  )
) +
  geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.20
  ) +
  scale_fill_manual(
    values = celltype_colors,
    limits = display_celltype_levels,
    breaks = display_celltype_levels,
    drop = FALSE
  ) +
  scale_x_discrete(labels = region_short_labels) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    tag = "C",
    title = NULL,
    x = NULL,
    y = "Proportion",
    fill = "Cell type"
  ) +
  theme_manuscript +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      face = "bold",
      size = 9
    )
  )

p1D <- ggplot(
  comparison_table[
    as.character(Final_celltype) %in%
      display_celltype_levels
  ],
  aes(
    x = Final_celltype,
    y = log2_response_nonresponse,
    fill = enriched_region
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey35",
    linewidth = 0.55
  ) +
  geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.20
  ) +
  scale_fill_manual(
    values = region_colors,
    labels = region_short_labels,
    drop = FALSE
  ) +
  labs(
    tag = "D",
    title = NULL,
    x = NULL,
    y = "log2 fold change\n(Response / Non-response)",
    fill = "Enriched region"
  ) +
  theme_manuscript +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 7.5
    )
  )

## ------------------------------------------------------------
## 8. Assemble and save Figure 1
## ------------------------------------------------------------

right_column <- (p1B / p1C / p1D) +
  plot_layout(
    heights = c(0.88, 1.08, 1.00),
    guides = "keep"
  )

figure1 <- (p1A | right_column) +
  plot_layout(
    widths = c(1.92, 1.08),
    guides = "keep"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0.005, 0.995)
  )

print(figure1)

save_plot_formats(
  plot_object = figure1,
  output_directory = figure1_dir,
  base_filename = "Figure1",
  width = 20.0,
  height = 11.4,
  save_tiff = TRUE
)

save_panel_collection(
  panel_list = list(p1A, p1B, p1C, p1D),
  output_directory = figure1_dir,
  prefix = "Figure1",
  width = c(13.0, 9.0, 8.0, 9.0),
  height = c(11.0, 6.5, 7.0, 6.5)
)

## ------------------------------------------------------------
## 9. Supplementary Figure 1: annotation marker DotPlot
## ------------------------------------------------------------

annotation_marker_sets <- list(
  Cancer = c(
    "EPCAM", "KRT8", "KRT18", "KRT19", "CEACAM6"
  ),
  CAF_apCAF = c(
    "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "CXCL12"
  ),
  CAF_ECM = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "FN1", "THBS2"
  ),
  CAF_iCAF = c(
    "IL6", "CXCL12", "CCL2", "IL1R1", "CXCL14"
  ),
  CAF_myCAF = c(
    "ACTA2", "TAGLN", "MYL9", "COL4A1", "COL4A2"
  ),
  Imm_Tcell = c(
    "PTPRC", "CD2", "CD3D", "CD3E", "IL7R"
  ),
  `Imm_Cytotoxic_NK/T` = c(
    "CD8A", "CD8B", "CCL5", "GZMH", "TIGIT"
  ),
  Imm_NK = c(
    "NKG7", "GNLY", "PRF1", "NCR1", "NCAM1"
  ),
  Imm_Myeloid = c(
    "CD68", "CD163", "CSF1R", "C1QA", "C1QB", "C1QC"
  ),
  Imm_Plasma = c(
    "JCHAIN", "MZB1", "IGKC", "IGHG1", "IGHA1"
  ),
  Endothelial = c(
    "PECAM1", "VWF", "EMCN", "KDR", "ENG"
  ),
  Pericyte = c(
    "RGS5", "PDGFRB", "CSPG4", "MCAM", "RBP1"
  )
)

available_features <- rownames(full_object)
used_markers <- character()

annotation_markers_available <- lapply(
  annotation_marker_sets,
  function(current_markers) {
    current_markers <- intersect(
      current_markers,
      available_features
    )

    current_markers <- current_markers[
      !current_markers %in% used_markers
    ]

    used_markers <<- c(
      used_markers,
      current_markers
    )

    current_markers
  }
)

annotation_markers_available <- annotation_markers_available[
  lengths(annotation_markers_available) > 0L
]

annotation_marker_table <- rbindlist(
  lapply(
    names(annotation_markers_available),
    function(current_celltype) {
      data.table(
        cell_type = current_celltype,
        marker = annotation_markers_available[[current_celltype]]
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  annotation_marker_table,
  file.path(
    supplementary_table1_dir,
    "Supplementary Table1_annotation_markers.csv"
  )
)

annotation_cells <- colnames(full_object)[
  as.character(full_object$Final_celltype) %in%
    celltype_levels
]

annotation_dot_object <- subset(
  full_object,
  cells = annotation_cells
)

annotation_dot_object$Final_celltype <- factor(
  as.character(annotation_dot_object$Final_celltype),
  levels = celltype_levels
)

Idents(annotation_dot_object) <- "Final_celltype"

pS1A <- DotPlot(
  annotation_dot_object,
  features = annotation_markers_available,
  assay = "Nanostring",
  group.by = "Final_celltype",
  cols = c("grey92", "#2166AC", "#B2182B"),
  dot.scale = 8,
  scale = TRUE,
  scale.by = "radius"
) +
  RotatedAxis() +
  labs(
    title = "Marker expression supporting final cell-type annotations",
    x = NULL,
    y = "Final cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      face = "bold",
      colour = "black",
      size = 8.2
    ),
    axis.text.y = element_text(
      face = "bold",
      colour = "black",
      size = 9.2
    ),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.22
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.55
    ),
    plot.margin = margin(8, 10, 8, 10)
  )

supplementary_figure1 <- pS1A +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0.005, 0.995)
  )

print(supplementary_figure1)

save_plot_formats(
  plot_object = supplementary_figure1,
  output_directory = supplementary1_dir,
  base_filename = "Supplementary Figure1",
  width = 18.0,
  height = 8.5,
  save_tiff = TRUE
)

save_plot_formats(
  plot_object = pS1A,
  output_directory = supplementary1_dir,
  base_filename = "Supplementary Figure1A",
  width = 18.0,
  height = 8.5,
  save_tiff = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    table1_dir,
    "Figure1_SupplementaryFigure1_sessionInfo.txt"
  )
)

cat("\n========================================\n")
cat("Figure 1 and Supplementary Figure 1 completed\n")
cat("Figure 1:\n", file.path(figure1_dir, "Figure1.pdf"), "\n")
cat(
  "Supplementary Figure 1:\n",
  file.path(
    supplementary1_dir,
    "Supplementary Figure1.pdf"
  ),
  "\n"
)
cat("========================================\n")

## ============================================================
## PART 2. Figure 2 and Supplementary Figure 2
## ============================================================

## ============================================================
## Figure 2 and Supplementary Figure 2
## Distance-dependent cancer-adjacent cell composition
##
## A-B Spatial maps of cancer and non-cancer cells within 25.39 um
## C   Cell-type composition across distance bands
## D   Response/non-response ratio within 25.39 um
## E   Distance-dependent proportions of CAF subtypes
## F   Distance-dependent proportions of immune subtypes
##
## Supplementary Figure 2
## A-B Aggregate CAF and immune proportions across radius thresholds
## C   Selected cell-type proportions across radius thresholds
## D   Response/non-response log2 ratio across radius thresholds
##
## Distance is measured from each non-cancer cell center to the
## nearest cancer-cell center in the same FOV. Each non-cancer cell
## is counted once. Cancer cells are excluded from every proportion.
## ============================================================

options(timeout = 1200, download.file.method = "libcurl")
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = 1200)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(FNN)
  library(scales)
})

set.seed(20261007)

## ------------------------------------------------------------
## 1. Paths and settings
## ------------------------------------------------------------

project_dir <- normalizePath(Sys.getenv("CRC_COSMX_PROJECT_DIR", unset = getwd()), mustWork = TRUE)

object_dir <- file.path(project_dir, "6 Analysis", "Data objects")
figure_dir <- file.path(project_dir, "6 Analysis", "Figures", "Figure2")
table_dir <- file.path(project_dir, "6 Analysis", "Tables", "Figure2")

supplementary2_dir <- file.path(
  project_dir,
  "6 Analysis",
  "Figures",
  "Supplementary Figure2"
)

supplementary2_table_dir <- file.path(
  project_dir,
  "6 Analysis",
  "Tables",
  "Supplementary Figure2"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(
  supplementary2_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  supplementary2_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

analysis_sample <- Sys.getenv("CRC_COSMX_SAMPLE", unset = "PUBLIC_SAMPLE_ID")
target_fovs <- c("4", "6")
pixel_size_um <- 0.12
primary_radius_um <- 25.39

region_mapping <- c(
  "4" = "Non-response region",
  "6" = "Response region"
)

region_levels <- c(
  "Non-response region",
  "Response region"
)

region_colors <- c(
  "Non-response region" = "#B2182B",
  "Response region" = "#2166AC"
)

distance_breaks <- c(
  -Inf,
  primary_radius_um,
  50,
  100,
  Inf
)

distance_labels <- c(
  paste0("≤", primary_radius_um, " µm"),
  paste0(">", primary_radius_um, "–50 µm"),
  ">50–100 µm",
  ">100 µm"
)

celltype_levels <- c(
  "CAF_apCAF",
  "CAF_ECM",
  "CAF_iCAF",
  "CAF_myCAF",
  "Imm_Tcell",
  "Imm_Cytotoxic_NK/T",
  "Imm_NK",
  "Imm_Myeloid",
  "Imm_Plasma",
  "Endothelial",
  "Pericyte"
)

caf_celltypes <- c(
  "CAF_apCAF",
  "CAF_ECM",
  "CAF_iCAF",
  "CAF_myCAF"
)

immune_celltypes <- c(
  "Imm_Tcell",
  "Imm_Cytotoxic_NK/T",
  "Imm_NK",
  "Imm_Myeloid",
  "Imm_Plasma"
)

## Final manuscript palette. Keep this identical in Figures 1-4.
celltype_colors <- c(
  "Cancer" = "#E3261C",
  "CAF_apCAF" = "#239B6E",
  "CAF_ECM" = "#6FC5AC",
  "CAF_iCAF" = "#53AC45",
  "CAF_myCAF" = "#28BF70",
  "Imm_Tcell" = "#C06BF0",
  "Imm_Cytotoxic_NK/T" = "#367EB5",
  "Imm_NK" = "#81BED8",
  "Imm_Myeloid" = "#5AA2CA",
  "Imm_Plasma" = "#ED79B7",
  "Endothelial" = "#686868",
  "Pericyte" = "#979797",
  "Unresolved" = "#BDBDBD"
)

celltype_colors <- celltype_colors[
  !duplicated(names(celltype_colors))
]

## Use the final annotated object already loaded in the session.
if (!exists("full_object")) {
  full_object_file <- file.path(
    object_dir,
    "CRC_CosMx_final_annotated_full_object.rds"
  )

  if (!file.exists(full_object_file)) {
    stop(
      paste0(
        "full_object is not loaded and the default RDS was not found: ",
        full_object_file,
        "\nLoad the final annotated object as full_object and rerun."
      )
    )
  }

  full_object <- readRDS(full_object_file)
}

required_metadata <- c(
  "sample_id",
  "fov",
  "Final_celltype",
  "CenterX_local_px",
  "CenterY_local_px"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(full_object[[]])
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

## ------------------------------------------------------------
## 2. Metadata and nearest-cancer distance
## ------------------------------------------------------------

spatial_data <- as.data.table(
  full_object[[]],
  keep.rownames = "cell"
)

spatial_data[, sample_id := as.character(sample_id)]
spatial_data[, fov := as.character(fov)]
spatial_data[, Final_celltype := as.character(Final_celltype)]
spatial_data[, x_px := as.numeric(CenterX_local_px)]
spatial_data[, y_px := as.numeric(CenterY_local_px)]

spatial_data <- spatial_data[
  sample_id == analysis_sample &
    fov %in% target_fovs &
    is.finite(x_px) &
    is.finite(y_px)
]

spatial_data[
  ,
  Region := factor(
    region_mapping[fov],
    levels = region_levels
  )
]

distance_result_list <- vector(
  mode = "list",
  length = length(target_fovs)
)

names(distance_result_list) <- target_fovs

for (current_fov in target_fovs) {
  current_data <- spatial_data[fov == current_fov]
  cancer_data <- current_data[Final_celltype == "Cancer"]
  noncancer_data <- current_data[Final_celltype != "Cancer"]

  if (nrow(cancer_data) == 0L) {
    stop("No cancer cells were found in FOV ", current_fov)
  }

  if (nrow(noncancer_data) == 0L) {
    stop("No non-cancer cells were found in FOV ", current_fov)
  }

  nearest_cancer <- FNN::get.knnx(
    data = as.matrix(
      cancer_data[, .(x_px, y_px)]
    ),
    query = as.matrix(
      noncancer_data[, .(x_px, y_px)]
    ),
    k = 1,
    algorithm = "kd_tree"
  )

  noncancer_data[
    ,
    nearest_cancer_index :=
      as.integer(nearest_cancer$nn.index[, 1])
  ]

  noncancer_data[
    ,
    nearest_cancer_cell :=
      cancer_data$cell[nearest_cancer_index]
  ]

  noncancer_data[
    ,
    distance_to_nearest_cancer_um :=
      as.numeric(nearest_cancer$nn.dist[, 1]) *
      pixel_size_um
  ]

  distance_result_list[[current_fov]] <- noncancer_data
}

noncancer_distance_data <- unique(
  rbindlist(
    distance_result_list,
    use.names = TRUE,
    fill = TRUE
  ),
  by = "cell"
)

if (
  any(
    noncancer_distance_data$Final_celltype == "Cancer",
    na.rm = TRUE
  )
) {
  stop("Cancer cells were detected in noncancer_distance_data.")
}

noncancer_distance_data[
  ,
  Distance_bin := cut(
    distance_to_nearest_cancer_um,
    breaks = distance_breaks,
    labels = distance_labels,
    include.lowest = TRUE,
    right = TRUE
  )
]

noncancer_distance_data[
  ,
  Distance_bin := factor(
    as.character(Distance_bin),
    levels = distance_labels,
    ordered = TRUE
  )
]

noncancer_distance_data[
  ,
  Final_celltype := factor(
    Final_celltype,
    levels = celltype_levels
  )
]

noncancer_distance_data <- noncancer_distance_data[
  !is.na(Final_celltype) &
    !is.na(Distance_bin)
]

cat("\nNon-cancer cells by region and distance:\n")

print(
  table(
    noncancer_distance_data$Region,
    noncancer_distance_data$Distance_bin
  )
)

## ------------------------------------------------------------
## 3. Cell-type proportions by distance
## ------------------------------------------------------------

composition_counts <- noncancer_distance_data[
  ,
  .(
    cells = .N
  ),
  by = .(
    Region,
    Distance_bin,
    Final_celltype
  )
]

composition_grid <- CJ(
  Region = region_levels,
  Distance_bin = distance_labels,
  Final_celltype = celltype_levels,
  unique = TRUE
)

composition_grid[
  ,
  Region := factor(
    Region,
    levels = region_levels
  )
]

composition_grid[
  ,
  Distance_bin := factor(
    Distance_bin,
    levels = distance_labels,
    ordered = TRUE
  )
]

composition_grid[
  ,
  Final_celltype := factor(
    Final_celltype,
    levels = celltype_levels
  )
]

composition_table <- merge(
  composition_grid,
  composition_counts,
  by = c(
    "Region",
    "Distance_bin",
    "Final_celltype"
  ),
  all.x = TRUE,
  sort = FALSE
)

composition_table[
  is.na(cells),
  cells := 0L
]

composition_table[
  ,
  total_noncancer_cells := sum(cells),
  by = .(
    Region,
    Distance_bin
  )
]

composition_table[
  ,
  proportion := fifelse(
    total_noncancer_cells > 0,
    cells / total_noncancer_cells,
    NA_real_
  )
]

## Jeffreys-style smoothing is used only for finite region ratios.
number_of_celltypes <- length(celltype_levels)

composition_table[
  ,
  smoothed_proportion :=
    (cells + 0.5) /
    (total_noncancer_cells + 0.5 * number_of_celltypes)
]

## ------------------------------------------------------------
## 4. Response/non-response comparison at every distance
## ------------------------------------------------------------

distance_comparison_wide <- dcast(
  composition_table,
  Distance_bin + Final_celltype ~ Region,
  value.var = c(
    "cells",
    "total_noncancer_cells",
    "proportion",
    "smoothed_proportion"
  )
)

distance_comparison_wide[
  ,
  log2_response_nonresponse_ratio := log2(
    `smoothed_proportion_Response region` /
      `smoothed_proportion_Non-response region`
  )
]

distance_comparison_wide[
  ,
  direction := factor(
    fifelse(
      log2_response_nonresponse_ratio >= 0,
      "Response region",
      "Non-response region"
    ),
    levels = region_levels
  )
]

## Primary cancer-adjacent band used in panel D.
near_label <- distance_labels[1]

near_wide <- copy(
  distance_comparison_wide[
    as.character(Distance_bin) == near_label
  ]
)

## ------------------------------------------------------------
## 5. Export tables
## ------------------------------------------------------------

fwrite(
  noncancer_distance_data[
    ,
    .(
      cell,
      sample_id,
      fov,
      Region,
      Final_celltype,
      nearest_cancer_cell,
      distance_to_nearest_cancer_um,
      Distance_bin
    )
  ],
  file.path(
    table_dir,
    "Figure2_cell_to_nearest_cancer_distance.csv"
  )
)

fwrite(
  composition_table,
  file.path(
    table_dir,
    "Figure2_celltype_proportion_by_distance.csv"
  )
)

fwrite(
  distance_comparison_wide,
  file.path(
    table_dir,
    "Figure2_distance_region_comparison.csv"
  )
)

fwrite(
  near_wide,
  file.path(
    table_dir,
    "Figure2_near_25um_region_comparison.csv"
  )
)

## ------------------------------------------------------------
## 6. Plot data for spatial panels
## ------------------------------------------------------------

plot_spatial_data <- copy(spatial_data)

plot_spatial_data[
  ,
  Plot_group := "Other non-cancer"
]

plot_spatial_data[
  Final_celltype == "Cancer",
  Plot_group := "Cancer"
]

near_celltype_lookup <- noncancer_distance_data[
  distance_to_nearest_cancer_um <= primary_radius_um,
  .(
    cell,
    near_celltype = as.character(Final_celltype)
  )
]

plot_spatial_data <- merge(
  plot_spatial_data,
  near_celltype_lookup,
  by = "cell",
  all.x = TRUE,
  sort = FALSE
)

plot_spatial_data[
  !is.na(near_celltype),
  Plot_group := near_celltype
]

plot_group_levels <- c(
  "Cancer",
  celltype_levels,
  "Other non-cancer"
)

plot_group_colors <- c(
  celltype_colors,
  "Other non-cancer" = "grey88"
)

plot_spatial_data[
  ,
  Plot_group := factor(
    Plot_group,
    levels = plot_group_levels
  )
]

## Draw background first, followed by nearby cells and cancer cells.
plot_spatial_data[
  ,
  plot_order := fifelse(
    Plot_group == "Other non-cancer",
    0L,
    fifelse(
      Plot_group == "Cancer",
      2L,
      1L
    )
  )
]

setorder(
  plot_spatial_data,
  plot_order
)

## ------------------------------------------------------------
## 7. Plotting helpers
## ------------------------------------------------------------

theme_figure2 <- theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    axis.title = element_text(
      colour = "black",
      face = "bold",
      size = 11
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    strip.text = element_text(
      face = "bold",
      colour = "black",
      size = 10
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(7, 9, 7, 9)
  )

make_spatial_plot <- function(region_name) {
  current_plot_data <- plot_spatial_data[
    Region == region_name
  ]

  ggplot(
    current_plot_data,
    aes(
      x = x_px * pixel_size_um,
      y = y_px * pixel_size_um,
      colour = Plot_group
    )
  ) +
    geom_point(
      size = 0.90,
      alpha = 0.88
    ) +
    scale_colour_manual(
      values = plot_group_colors,
      drop = FALSE
    ) +
    coord_equal() +
    labs(
      title = region_name,
      subtitle = paste0(
        "Colored non-cancer cells: ≤",
        primary_radius_um,
        " µm from cancer"
      ),
      x = "X (µm)",
      y = "Y (µm)",
      colour = "Cell type"
    ) +
    theme_figure2 +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 10
      )
    )
}

make_proportion_trend_plot <- function(
  selected_celltypes,
  facet_columns
) {
  trend_data <- copy(
    composition_table[
    as.character(Final_celltype) %in%
      selected_celltypes
    ]
  )

  trend_data[
    ,
    Final_celltype := factor(
      as.character(Final_celltype),
      levels = selected_celltypes
    )
  ]

  ggplot(
    trend_data,
    aes(
      x = Distance_bin,
      y = proportion,
      colour = Region,
      group = Region
    )
  ) +
    geom_line(
      linewidth = 1.05,
      na.rm = TRUE
    ) +
    geom_point(
      size = 2.7,
      stroke = 0.25,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = region_colors,
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 0.1
      ),
      expand = expansion(
        mult = c(0.05, 0.12)
      )
    ) +
    facet_wrap(
      vars(Final_celltype),
      ncol = facet_columns,
      scales = "free_y"
    ) +
    labs(
      x = "Distance to nearest cancer cell",
      y = "Cell-type proportion",
      colour = "Region"
    ) +
    theme_figure2 +
    theme(
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        face = "bold",
        size = 8
      ),
      strip.text = element_text(
        face = "bold",
        size = 8.5
      ),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(
        size = 8
      ),
      panel.grid.major.y = element_line(
        colour = "grey90",
        linewidth = 0.25
      )
    )
}

## ------------------------------------------------------------
## 8. Build panels A-F
## ------------------------------------------------------------

pA <- make_spatial_plot(
  "Non-response region"
)

pB <- make_spatial_plot(
  "Response region"
)

pC <- ggplot(
  composition_table,
  aes(
    x = Distance_bin,
    y = proportion,
    fill = Final_celltype
  )
) +
  geom_col(
    width = 0.78,
    colour = "white",
    linewidth = 0.15
  ) +
  scale_fill_manual(
    values = celltype_colors,
    limits = c(
      "Cancer",
      celltype_levels
    ),
    breaks = c(
      "Cancer",
      celltype_levels
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  facet_wrap(
    vars(Region),
    ncol = 1
  ) +
  labs(
    title = "Non-cancer cell composition by distance",
    x = "Distance to nearest cancer cell",
    y = "Cell-type proportion",
    fill = "Cell type"
  ) +
  theme_figure2 +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold",
      size = 8
    ),
    legend.position = "right",
    legend.text = element_text(
      size = 8
    )
  )

near_wide[
  ,
  Final_celltype := factor(
    as.character(Final_celltype),
    levels = rev(celltype_levels)
  )
]

pD <- ggplot(
  near_wide,
  aes(
    x = log2_response_nonresponse_ratio,
    y = Final_celltype,
    colour = direction
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.6,
    colour = "grey35"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2_response_nonresponse_ratio,
      yend = Final_celltype
    ),
    linewidth = 0.85
  ) +
  geom_point(
    size = 3.2
  ) +
  scale_colour_manual(
    values = region_colors,
    drop = FALSE
  ) +
  labs(
    x = "log2(Response / Non-response proportion)",
    y = NULL,
    colour = "Higher in"
  ) +
  theme_figure2 +
  theme(
    axis.text.y = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "top"
  )

pE <- make_proportion_trend_plot(
  selected_celltypes = caf_celltypes,
  facet_columns = 2
)

pF <- make_proportion_trend_plot(
  selected_celltypes = immune_celltypes,
  facet_columns = 2
)

figure2 <- (pA | pB | pC) /
  (pD | pE | pF) +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(
    widths = c(1, 1, 1.18),
    heights = c(1, 1),
    guides = "keep"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0.01, 0.99)
  )

print(figure2)

## ------------------------------------------------------------
## 9. Save complete figure and individual panels
## ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2.pdf"
  ),
  plot = figure2,
  width = 21,
  height = 12,
  units = "in",
  device = cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2.png"
  ),
  plot = figure2,
  width = 21,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2.tiff"
  ),
  plot = figure2,
  width = 21,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

panel_list <- list(
  Figure2A = pA,
  Figure2B = pB,
  Figure2C = pC,
  Figure2D = pD,
  Figure2E = pE,
  Figure2F = pF
)

for (panel_name in names(panel_list)) {
  panel_plot <- panel_list[[panel_name]]

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(panel_name, ".pdf")
    ),
    plot = panel_plot,
    width = 8.8,
    height = 6.8,
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(panel_name, ".png")
    ),
    plot = panel_plot,
    width = 8.8,
    height = 6.8,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

cat("\n========================================\n")
cat("Figure 2 distance-composition analysis completed\n")
cat("Primary cancer-adjacent radius:", primary_radius_um, "µm\n")
cat("\nNear-region comparison:\n")

print(
  near_wide[
    ,
    .(
      Final_celltype,
      `Non-response cells` =
        `cells_Non-response region`,
      `Response cells` =
        `cells_Response region`,
      log2_ratio =
        log2_response_nonresponse_ratio
    )
  ]
)

cat("\nFigures:", figure_dir, "\n")
cat("Tables:", table_dir, "\n")

## ------------------------------------------------------------
## 10. Supplementary Figure 2: distance-threshold sensitivity
## ------------------------------------------------------------

sensitivity_radii_um <- c(
  12.69,
  25.39,
  38.08,
  50.77,
  101.56
)

sensitivity_list <- vector(
  mode = "list",
  length = length(sensitivity_radii_um)
)

for (radius_index in seq_along(sensitivity_radii_um)) {
  current_radius <- sensitivity_radii_um[radius_index]

  current_counts <- noncancer_distance_data[
    distance_to_nearest_cancer_um <= current_radius,
    .(cell_count = .N),
    by = .(
      Region,
      Final_celltype
    )
  ]

  current_complete <- merge(
    CJ(
      Region = factor(region_levels, levels = region_levels),
      Final_celltype = factor(
        celltype_levels,
        levels = celltype_levels
      ),
      unique = TRUE
    ),
    current_counts,
    by = c("Region", "Final_celltype"),
    all.x = TRUE,
    sort = FALSE
  )

  current_complete[
    is.na(cell_count),
    cell_count := 0L
  ]

  current_complete[
    ,
    total_noncancer_cells := sum(cell_count),
    by = Region
  ]

  current_complete[
    ,
    proportion := fifelse(
      total_noncancer_cells > 0L,
      cell_count / total_noncancer_cells,
      NA_real_
    )
  ]

  current_complete[, radius_um := current_radius]

  sensitivity_list[[radius_index]] <- current_complete
}

sensitivity_data <- rbindlist(
  sensitivity_list,
  use.names = TRUE,
  fill = TRUE
)

sensitivity_data[
  ,
  radius_label := factor(
    paste0(
      "≤",
      formatC(
        radius_um,
        format = "f",
        digits = 2
      ),
      " µm"
    ),
    levels = paste0(
      "≤",
      formatC(
        sensitivity_radii_um,
        format = "f",
        digits = 2
      ),
      " µm"
    )
  )
]

aggregate_sensitivity <- sensitivity_data[
  ,
  .(
    proportion = sum(proportion, na.rm = TRUE)
  ),
  by = .(
    Region,
    radius_um,
    radius_label,
    compartment = fifelse(
      as.character(Final_celltype) %in% caf_celltypes,
      "CAF",
      fifelse(
        as.character(Final_celltype) %in% immune_celltypes,
        "Immune",
        "Other"
      )
    )
  )
][
  compartment %in% c("CAF", "Immune")
]

sensitivity_wide <- dcast(
  sensitivity_data,
  radius_um + radius_label + Final_celltype ~ Region,
  value.var = "proportion",
  fill = 0
)

sensitivity_wide[
  ,
  log2_response_nonresponse_ratio := log2(
    (`Response region` + 1e-6) /
      (`Non-response region` + 1e-6)
  )
]

sensitivity_wide[
  ,
  Final_celltype := factor(
    Final_celltype,
    levels = rev(celltype_levels)
  )
]

fwrite(
  sensitivity_data,
  file.path(
    supplementary2_table_dir,
    "SupplementaryFigure2_celltype_proportion_by_radius.csv"
  )
)

fwrite(
  sensitivity_wide,
  file.path(
    supplementary2_table_dir,
    "SupplementaryFigure2_response_nonresponse_ratio_by_radius.csv"
  )
)

make_aggregate_sensitivity_plot <- function(
  compartment_name,
  plot_title
) {
  current_data <- aggregate_sensitivity[
    compartment == compartment_name
  ]

  ggplot(
    current_data,
    aes(
      x = radius_um,
      y = proportion,
      colour = Region,
      group = Region
    )
  ) +
    geom_line(linewidth = 1.15) +
    geom_point(size = 3.0) +
    scale_colour_manual(
      values = region_colors,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = sensitivity_radii_um,
      labels = formatC(
        sensitivity_radii_um,
        format = "f",
        digits = 2
      )
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 0.1),
      expand = expansion(mult = c(0.04, 0.12))
    ) +
    labs(
      title = plot_title,
      x = "Maximum distance to nearest cancer cell (µm)",
      y = "Proportion among non-cancer cells",
      colour = NULL
    ) +
    theme_figure2 +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_line(
        colour = "grey90",
        linewidth = 0.25
      )
    )
}

pS2A <- make_aggregate_sensitivity_plot(
  compartment_name = "CAF",
  plot_title = "CAF proportion across distance thresholds"
)

pS2B <- make_aggregate_sensitivity_plot(
  compartment_name = "Immune",
  plot_title = "Immune proportion across distance thresholds"
)

selected_sensitivity_celltypes <- c(
  "CAF_ECM",
  "Imm_Tcell",
  "Imm_Cytotoxic_NK/T",
  "Imm_Myeloid",
  "Imm_Plasma"
)

pS2C_data <- copy(
  sensitivity_data[
    as.character(Final_celltype) %in%
      selected_sensitivity_celltypes
  ]
)

pS2C_data[
  ,
  Final_celltype := factor(
    as.character(Final_celltype),
    levels = selected_sensitivity_celltypes
  )
]

pS2C <- ggplot(
  pS2C_data,
  aes(
    x = radius_um,
    y = proportion,
    colour = Region,
    group = Region
  )
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.4) +
  scale_colour_manual(
    values = region_colors,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = sensitivity_radii_um,
    labels = formatC(
      sensitivity_radii_um,
      format = "f",
      digits = 2
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  facet_wrap(
    vars(Final_celltype),
    ncol = 2,
    scales = "free_y"
  ) +
  labs(
    title = "Selected cell populations across distance thresholds",
    x = "Maximum distance to nearest cancer cell (µm)",
    y = "Cell-type proportion",
    colour = NULL
  ) +
  theme_figure2 +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 8.5),
    panel.grid.major.y = element_line(
      colour = "grey90",
      linewidth = 0.25
    )
  )

pS2D <- ggplot(
  sensitivity_wide,
  aes(
    x = radius_label,
    y = Final_celltype,
    fill = log2_response_nonresponse_ratio
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.55
  ) +
  geom_text(
    aes(
      label = formatC(
        log2_response_nonresponse_ratio,
        format = "f",
        digits = 2
      )
    ),
    size = 2.7,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = region_colors[["Non-response region"]],
    mid = "white",
    high = region_colors[["Response region"]],
    midpoint = 0,
    oob = squish,
    name = "log2(Response /\nNon-response)"
  ) +
  labs(
    title = "Regional composition difference across distance thresholds",
    x = "Maximum distance to nearest cancer cell",
    y = NULL
  ) +
  theme_figure2 +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 8
    ),
    legend.position = "right",
    panel.grid = element_blank()
  )

supplementary_figure2 <- wrap_plots(
  pS2A,
  pS2B,
  pS2C,
  pS2D,
  nrow = 1,
  widths = c(1, 1, 1.30, 1.20),
  guides = "keep"
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0.005, 0.995)
  )

print(supplementary_figure2)

ggsave(
  filename = file.path(
    supplementary2_dir,
    "Supplementary Figure2.pdf"
  ),
  plot = supplementary_figure2,
  width = 25,
  height = 7.5,
  units = "in",
  device = cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary2_dir,
    "Supplementary Figure2.png"
  ),
  plot = supplementary_figure2,
  width = 25,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary2_dir,
    "Supplementary Figure2.tiff"
  ),
  plot = supplementary_figure2,
  width = 25,
  height = 7.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

supplementary2_panels <- list(
  pS2A,
  pS2B,
  pS2C,
  pS2D
)

supplementary2_panel_widths <- c(8.0, 8.0, 11.0, 10.0)
supplementary2_panel_heights <- c(6.2, 6.2, 8.0, 7.5)

for (panel_index in seq_along(supplementary2_panels)) {
  panel_letter <- LETTERS[panel_index]

  ggsave(
    filename = file.path(
      supplementary2_dir,
      paste0("Supplementary Figure2", panel_letter, ".pdf")
    ),
    plot = supplementary2_panels[[panel_index]],
    width = supplementary2_panel_widths[panel_index],
    height = supplementary2_panel_heights[panel_index],
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(
      supplementary2_dir,
      paste0("Supplementary Figure2", panel_letter, ".png")
    ),
    plot = supplementary2_panels[[panel_index]],
    width = supplementary2_panel_widths[panel_index],
    height = supplementary2_panel_heights[panel_index],
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(
    table_dir,
    "Figure2_SupplementaryFigure2_sessionInfo.txt"
  )
)

cat(
  "Supplementary Figure 2:",
  file.path(
    supplementary2_dir,
    "Supplementary Figure2.pdf"
  ),
  "\n"
)
cat("========================================\n")

## ============================================================
## PART 3. Figure 3, Figure 4, and Supplementary Figures 3-4
## ============================================================

## ============================================================
## Figure 3 and Figure 4: post-radiotherapy regional ecosystems
##
## Study design
## - One patient with locally advanced rectal adenocarcinoma
## - Post-radiotherapy surgical specimen
## - Non-response region: FOV4
## - Response region: FOV6
## - Primary tumor-adjacent radius: 25.39 um
##
## Figure 3
## A-B  Cancer, ECM-CAF, and myeloid maps within 101.56 um
## C    ECM-CAF matrix-remodeling score comparison
## D    Myeloid inflammatory score comparison
## E    ECM-CAF abundance across distance from cancer cells
## F    ECM-CAF remodeling score across continuous cancer distance
## G    Myeloid-cell abundance across distance from cancer cells
## H    Myeloid inflammatory score across continuous cancer distance
## I-J  Cell-type-specific gene-expression effect sizes
## K    Coordinated CAF-myeloid distance heatmap
##
## Supplementary Figure 3
## A-B  T-cell and cytotoxic NK/T functional-state maps
## C-D  T-cell and cytotoxic NK/T score comparisons
## E-F  T-cell and cytotoxic NK/T gene-expression effect sizes
##
## Figure 4
## A-B  Integrated cancer IFN, ECM-CAF, and myeloid state maps
## C    Cancer IFN/adaptive-resistance score comparison
## D    Cancer antigen-presentation score comparison
## E    Selected cancer-cell gene effect sizes
## F    Selected cancer-centered spatial ligand-receptor co-detection
##
## Supplementary Figure 4
## A    Extended cancer-centered ligand-receptor co-detection
## B    Cancer IFN/adaptive resistance versus local ECM-CAF burden
##
## Important interpretation
## - This is a within-patient, one-FOV-per-region comparison.
## - Cell-level raw P values are exploratory and are not biological-replicate P values.
## - The plots describe post-RT regional associations, not treatment-induced causality.
## ============================================================

options(timeout = 1200, download.file.method = "libcurl")
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = 1200)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(FNN)
  library(Matrix)
  library(scales)
})

set.seed(20261008)

## ------------------------------------------------------------
## 1. Paths and settings
## ------------------------------------------------------------

project_dir <- normalizePath(Sys.getenv("CRC_COSMX_PROJECT_DIR", unset = getwd()), mustWork = TRUE)

analysis_dir <- file.path(project_dir, "6 Analysis")
object_dir <- file.path(analysis_dir, "Data objects")

figure3_dir <- file.path(analysis_dir, "Figures", "Figure3")
figure4_dir <- file.path(analysis_dir, "Figures", "Figure4")
supplementary3_dir <- file.path(
  analysis_dir,
  "Figures",
  "Supplementary Figure3"
)
supplementary4_dir <- file.path(
  analysis_dir,
  "Figures",
  "Supplementary Figure4"
)

table3_dir <- file.path(analysis_dir, "Tables", "Figure3")
table4_dir <- file.path(analysis_dir, "Tables", "Figure4")

output_directories <- c(
  figure3_dir,
  figure4_dir,
  supplementary3_dir,
  supplementary4_dir,
  table3_dir,
  table4_dir
)

for (current_directory in output_directories) {
  dir.create(
    current_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

## Set this manually only when the annotated object has a different name.
full_object_file <- file.path(
  object_dir,
  "CRC_CosMx_final_annotated_full_object.rds"
)

analysis_sample <- Sys.getenv("CRC_COSMX_SAMPLE", unset = "PUBLIC_SAMPLE_ID")
target_fovs <- c("4", "6")
pixel_size_um <- 0.12
primary_radius_um <- 25.39
niche_radius_um <- primary_radius_um * 4
minimum_interaction_sender_cells <- 10L
minimum_interaction_receiver_cells <- 3L
minimum_codetection_fraction <- 0.01
maximum_main_interactions <- 12L
interaction_pseudocount <- 1e-6

region_mapping <- c(
  "4" = "Non-response region",
  "6" = "Response region"
)

region_levels <- c(
  "Non-response region",
  "Response region"
)

region_colors <- c(
  "Non-response region" = "#B2182B",
  "Response region" = "#2166AC"
)

## Final manuscript palette. Keep this identical in Figures 1-4.
celltype_colors <- c(
  "Cancer" = "#E3261C",
  "CAF_apCAF" = "#239B6E",
  "CAF_ECM" = "#6FC5AC",
  "CAF_iCAF" = "#53AC45",
  "CAF_myCAF" = "#28BF70",
  "Imm_Tcell" = "#C06BF0",
  "Imm_Cytotoxic_NK/T" = "#367EB5",
  "Imm_NK" = "#81BED8",
  "Imm_Myeloid" = "#5AA2CA",
  "Imm_Plasma" = "#ED79B7",
  "Endothelial" = "#686868",
  "Pericyte" = "#979797",
  "Unresolved" = "#BDBDBD"
)

celltype_colors <- celltype_colors[
  !duplicated(names(celltype_colors))
]

score_low_color <- "#2166AC"
score_high_color <- "#B2182B"

minimum_markers_per_program <- 3L
module_score_ctrl <- 25L
module_score_nbin <- 12L
save_scored_object <- TRUE

## ------------------------------------------------------------
## 2. Marker programs used in the main and supplementary figures
## ------------------------------------------------------------

marker_sets <- list(
  CAF_ECM_remodeling = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2",
    "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL11A1",
    "COL12A1", "FN1", "DCN", "LUM", "THBS2", "POSTN",
    "SPARC", "MMP1", "MMP2", "MMP11", "MRC2", "VCAN",
    "FAP", "INHBA"
  ),

  Tcell_effector_activation = c(
    "CD69", "ICOS", "TNFRSF9", "IL2RA", "IFNG", "GZMK",
    "CCL5", "IL32", "RUNX3", "TBX21"
  ),

  Cytotoxic_NKT_cytotoxicity = c(
    "CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "GZMH",
    "GZMK", "CCL5", "CTSW", "KLRD1", "KLRK1", "KLRF1",
    "FGFBP2", "EOMES", "TBX21", "IFNG"
  ),

  Cytotoxic_NKT_exhaustion = c(
    "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4",
    "ENTPD1"
  ),

  Myeloid_inflammatory = c(
    "IL1B", "CXCL8", "S100A8", "S100A9", "TNF", "NFKBIA",
    "IL1RN", "CCL3/L1/L3", "PTGS2", "HCAR2/3", "CSF3R"
  ),

  Myeloid_antigen_presentation = c(
    "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
    "HLA-DQB1/2", "CD74", "CIITA"
  ),

  Plasma_antibody_secretion = c(
    "JCHAIN", "MZB1", "IGKC", "IGHG1", "IGHG2", "IGHA1",
    "SDC1", "XBP1", "PRDM1"
  ),

  Cancer_checkpoint_ligand = c(
    "CD274", "PDCD1LG2", "PVR", "NECTIN2", "LGALS9", "CD47",
    "TNFRSF14", "HLA-E"
  ),

  Cancer_antigen_presentation = c(
    "B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "NLRC5",
    "PSMB8", "PSMB9"
  ),

  Cancer_immunosuppressive = c(
    "TGFB1", "IDO1", "NT5E", "ENTPD1", "CSF1", "VEGFA", "IL10"
  ),

  Cancer_IFN_adaptive_resistance = c(
    "STAT1", "IRF1", "ISG15", "IFIT1", "IFIT3", "MX1", "CXCL9",
    "CXCL10", "B2M", "CD274", "IDO1"
  )
)

functional_definitions <- data.table(
  program = c(
    "CAF_ECM_remodeling",
    "Tcell_effector_activation",
    "Cytotoxic_NKT_cytotoxicity",
    "Cytotoxic_NKT_exhaustion",
    "Myeloid_inflammatory",
    "Myeloid_antigen_presentation",
    "Plasma_antibody_secretion"
  ),
  target_celltype = c(
    "CAF_ECM",
    "Imm_Tcell",
    "Imm_Cytotoxic_NK/T",
    "Imm_Cytotoxic_NK/T",
    "Imm_Myeloid",
    "Imm_Myeloid",
    "Imm_Plasma"
  ),
  display_label = c(
    "ECM-CAF matrix remodeling",
    "T-cell effector/activation",
    "Cytotoxic NK/T cytotoxicity",
    "Cytotoxic NK/T exhaustion",
    "Myeloid inflammatory activity",
    "Myeloid antigen presentation",
    "Plasma-cell antibody secretion"
  )
)

cancer_programs <- c(
  "Cancer_checkpoint_ligand",
  "Cancer_antigen_presentation",
  "Cancer_immunosuppressive",
  "Cancer_IFN_adaptive_resistance"
)

cancer_program_display <- c(
  "Cancer_checkpoint_ligand" = "Checkpoint ligands",
  "Cancer_antigen_presentation" = "Antigen presentation",
  "Cancer_immunosuppressive" = "Immunosuppressive signaling",
  "Cancer_IFN_adaptive_resistance" = "IFN/adaptive resistance"
)

selected_gene_groups <- list(
  "IRDS-like program" = c(
    "STAT1", "ISG15", "IFIT1", "IFIT3", "MX1"
  ),
  "Suppressive/myeloid recruitment" = c(
    "IDO1", "CSF1", "CD47", "TNFRSF14"
  ),
  "Antigen presentation" = c(
    "B2M", "NLRC5", "TAP1", "TAP2"
  ),
  "Immune recruitment" = c(
    "CXCL9", "CXCL10"
  )
)

## ------------------------------------------------------------
## 3. Reusable helper functions
## ------------------------------------------------------------

safe_wilcoxon <- function(x, y) {
  x <- as.numeric(x[is.finite(x)])
  y <- as.numeric(y[is.finite(y)])

  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }

  suppressWarnings(
    wilcox.test(
      x = x,
      y = y,
      alternative = "two.sided",
      exact = FALSE
    )$p.value
  )
}

format_raw_p <- function(value) {
  if (!is.finite(value)) {
    return("P = NA")
  }

  if (value < 0.001) {
    return("P < 0.001")
  }

  paste0("P = ", formatC(value, format = "f", digits = 3))
}

raw_p_symbol <- function(value) {
  if (!is.finite(value)) {
    return("")
  }

  if (value < 0.001) {
    return("***")
  }

  if (value < 0.01) {
    return("**")
  }

  if (value < 0.05) {
    return("*")
  }

  ""
}

safe_spearman <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  complete_rows <- is.finite(x) & is.finite(y)
  x <- x[complete_rows]
  y <- y[complete_rows]

  if (
    length(x) < 3L ||
      length(unique(x)) < 2L ||
      length(unique(y)) < 2L
  ) {
    return(
      list(
        n = as.integer(length(x)),
        rho = NA_real_,
        raw_p = NA_real_
      )
    )
  }

  test_result <- suppressWarnings(
    cor.test(
      x = x,
      y = y,
      method = "spearman",
      exact = FALSE
    )
  )

  list(
    n = as.integer(length(x)),
    rho = as.numeric(unname(test_result$estimate)),
    raw_p = as.numeric(test_result$p.value)
  )
}

standardize_vector <- function(values) {
  values <- as.numeric(values)
  current_sd <- sd(values, na.rm = TRUE)

  if (!is.finite(current_sd) || current_sd == 0) {
    return(rep(0, length(values)))
  }

  as.numeric(
    (values - mean(values, na.rm = TRUE)) / current_sd
  )
}

score_limit <- function(values) {
  values <- abs(as.numeric(values[is.finite(values)]))

  if (length(values) == 0L) {
    return(1)
  }

  current_limit <- as.numeric(
    quantile(
      values,
      probs = 0.98,
      na.rm = TRUE,
      names = FALSE
    )
  )

  if (!is.finite(current_limit) || current_limit <= 0) {
    current_limit <- 1
  }

  current_limit
}

safe_rank_biserial <- function(
  response_values,
  nonresponse_values
) {
  response_values <- as.numeric(
    response_values[is.finite(response_values)]
  )

  nonresponse_values <- as.numeric(
    nonresponse_values[is.finite(nonresponse_values)]
  )

  response_n <- length(response_values)
  nonresponse_n <- length(nonresponse_values)

  if (response_n < 2L || nonresponse_n < 2L) {
    return(NA_real_)
  }

  combined_values <- c(
    response_values,
    nonresponse_values
  )

  combined_ranks <- rank(
    combined_values,
    ties.method = "average"
  )

  response_rank_sum <- sum(
    combined_ranks[seq_len(response_n)]
  )

  mann_whitney_u <- response_rank_sum -
    response_n * (response_n + 1) / 2

  as.numeric(
    2 * mann_whitney_u /
      (response_n * nonresponse_n) - 1
  )
}

calculate_sender_local_interaction <- function(
  sender_table,
  receiver_table,
  expression_matrix,
  ligand_gene,
  receptor_gene,
  radius_um
) {
  sender_cells <- intersect(
    sender_table$cell,
    colnames(expression_matrix)
  )

  receiver_cells <- intersect(
    receiver_table$cell,
    colnames(expression_matrix)
  )

  output <- copy(
    sender_table[match(sender_cells, cell)]
  )

  output[, ligand_expression := 0]
  output[, local_receiver_count := 0L]
  output[, local_mean_receptor_expression := 0]
  output[, local_receptor_positive_fraction := 0]
  output[, local_interaction_score := 0]
  output[, ligand_receptor_codetected := FALSE]

  if (nrow(output) == 0L) {
    return(output)
  }

  ligand_expression <- as.numeric(
    expression_matrix[
      ligand_gene,
      sender_cells,
      drop = TRUE
    ]
  )

  output[, ligand_expression := ligand_expression]

  if (length(receiver_cells) == 0L) {
    return(output)
  }

  receiver_table <- copy(
    receiver_table[match(receiver_cells, cell)]
  )

  receptor_expression <- as.numeric(
    expression_matrix[
      receptor_gene,
      receiver_cells,
      drop = TRUE
    ]
  )

  nearest_result <- FNN::get.knnx(
    data = as.matrix(receiver_table[, .(x_um, y_um)]),
    query = as.matrix(output[, .(x_um, y_um)]),
    k = nrow(receiver_table),
    algorithm = "kd_tree"
  )

  neighbor_indices <- nearest_result$nn.index
  neighbor_distances <- nearest_result$nn.dist

  if (is.null(dim(neighbor_indices))) {
    neighbor_indices <- matrix(neighbor_indices, ncol = 1L)
    neighbor_distances <- matrix(neighbor_distances, ncol = 1L)
  }

  local_counts <- integer(nrow(output))
  local_receptor_means <- numeric(nrow(output))
  local_receptor_positive_fractions <- numeric(nrow(output))
  local_scores <- numeric(nrow(output))
  local_codetection <- logical(nrow(output))

  for (cell_index in seq_len(nrow(output))) {
    inside_radius <- which(
      neighbor_distances[cell_index, ] <= radius_um
    )

    local_counts[cell_index] <- length(inside_radius)

    if (length(inside_radius) > 0L) {
      receiver_rows <- neighbor_indices[
        cell_index,
        inside_radius,
        drop = TRUE
      ]

      local_receptor_values <- receptor_expression[
        receiver_rows
      ]

      local_receptor_means[cell_index] <- mean(
        local_receptor_values,
        na.rm = TRUE
      )

      local_receptor_positive_fractions[cell_index] <- mean(
        local_receptor_values > 0,
        na.rm = TRUE
      )

      local_scores[cell_index] <- sqrt(
        pmax(ligand_expression[cell_index], 0) *
          pmax(local_receptor_means[cell_index], 0)
      )

      local_codetection[cell_index] <-
        ligand_expression[cell_index] > 0 &&
        any(local_receptor_values > 0)
    }
  }

  output[, local_receiver_count := as.integer(local_counts)]
  output[
    ,
    local_mean_receptor_expression :=
      as.numeric(local_receptor_means)
  ]
  output[
    ,
    local_receptor_positive_fraction :=
      as.numeric(local_receptor_positive_fractions)
  ]
  output[, local_interaction_score := as.numeric(local_scores)]
  output[
    ,
    ligand_receptor_codetected := as.logical(local_codetection)
  ]

  output
}

theme_manuscript <- theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(
      colour = "black",
      face = "bold",
      size = 8.5
    ),
    axis.title = element_text(
      colour = "black",
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 11,
      margin = margin(b = 3)
    ),
    plot.subtitle = element_blank(),
    strip.text = element_text(
      face = "bold",
      colour = "black",
      size = 9
    ),
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(size = 8),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.55
    ),
    plot.margin = margin(6, 7, 6, 7)
  )

save_plot_formats <- function(
  plot_object,
  output_directory,
  base_filename,
  width,
  height,
  save_tiff = FALSE
) {
  ggsave(
    filename = file.path(output_directory, paste0(base_filename, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(output_directory, paste0(base_filename, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  if (isTRUE(save_tiff)) {
    ggsave(
      filename = file.path(output_directory, paste0(base_filename, ".tiff")),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
  }
}

save_panel_collection <- function(
  panel_list,
  output_directory,
  prefix,
  width,
  height
) {
  for (panel_index in seq_along(panel_list)) {
    panel_letter <- LETTERS[panel_index]

    save_plot_formats(
      plot_object = panel_list[[panel_index]],
      output_directory = output_directory,
      base_filename = paste0(prefix, panel_letter),
      width = width[panel_index],
      height = height[panel_index],
      save_tiff = FALSE
    )
  }
}

## ------------------------------------------------------------
## 4. Load and validate the final annotated object
## ------------------------------------------------------------

if (!exists("full_object")) {
  if (file.exists(full_object_file)) {
    full_object <- readRDS(full_object_file)
    message("Loaded annotated object: ", full_object_file)
  } else {
    candidate_files <- list.files(
      object_dir,
      pattern = "\\.rds$",
      full.names = TRUE
    )

    matched_file <- NA_character_

    for (candidate_file in rev(candidate_files)) {
      candidate_object <- tryCatch(
        readRDS(candidate_file),
        error = function(e) NULL
      )

      is_valid_object <- inherits(candidate_object, "Seurat") &&
        "Final_celltype" %in% colnames(candidate_object[[]])

      if (is_valid_object) {
        full_object <- candidate_object
        matched_file <- candidate_file
        break
      }

      rm(candidate_object)
      invisible(gc())
    }

    if (is.na(matched_file)) {
      stop(
        "No annotated Seurat RDS containing Final_celltype was found. ",
        "Set full_object_file manually."
      )
    }

    message("Loaded annotated object: ", matched_file)
  }
}

if (!inherits(full_object, "Seurat")) {
  stop("full_object is not a Seurat object.")
}

if (!"Nanostring" %in% Assays(full_object)) {
  stop("The Nanostring assay is missing from full_object.")
}

DefaultAssay(full_object) <- "Nanostring"

required_metadata <- c(
  "sample_id",
  "fov",
  "Final_celltype",
  "CenterX_local_px",
  "CenterY_local_px"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(full_object[[]])
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

## ------------------------------------------------------------
## 5. Calculate module scores once and retain stable column names
## ------------------------------------------------------------

available_features <- rownames(full_object[["Nanostring"]])
score_columns <- character(0)
marker_availability_list <- vector(
  mode = "list",
  length = length(marker_sets)
)
names(marker_availability_list) <- names(marker_sets)

for (current_program in names(marker_sets)) {
  requested_markers <- unique(marker_sets[[current_program]])
  current_markers <- intersect(
    requested_markers,
    available_features
  )

  marker_availability_list[[current_program]] <- data.table(
    program = current_program,
    gene = requested_markers,
    available = requested_markers %in% current_markers
  )

  if (length(current_markers) < minimum_markers_per_program) {
    warning(
      "Skipping ",
      current_program,
      ": only ",
      length(current_markers),
      " marker(s) are available."
    )
    next
  }

  safe_program <- gsub(
    pattern = "[^A-Za-z0-9]+",
    replacement = "_",
    x = current_program
  )

  temporary_prefix <- paste0("TmpPostRT_", safe_program, "_")
  final_column <- paste0("PostRT_", safe_program, "_score")

  if (final_column %in% colnames(full_object[[]])) {
    full_object[[final_column]] <- NULL
  }

  old_temporary_columns <- grep(
    pattern = paste0("^", temporary_prefix),
    x = colnames(full_object[[]]),
    value = TRUE
  )

  if (length(old_temporary_columns) > 0L) {
    for (old_column in old_temporary_columns) {
      full_object[[old_column]] <- NULL
    }
  }

  full_object <- AddModuleScore(
    object = full_object,
    features = list(current_markers),
    assay = "Nanostring",
    name = temporary_prefix,
    ctrl = module_score_ctrl,
    nbin = module_score_nbin,
    seed = 20261008,
    search = FALSE
  )

  generated_columns <- grep(
    pattern = paste0("^", temporary_prefix),
    x = colnames(full_object[[]]),
    value = TRUE
  )

  if (length(generated_columns) != 1L) {
    stop(
      "Unexpected score columns for ",
      current_program,
      ": ",
      paste(generated_columns, collapse = ", ")
    )
  }

  full_object[[final_column]] <- as.numeric(
    full_object[[generated_columns[1], drop = TRUE]]
  )

  full_object[[generated_columns[1]]] <- NULL
  score_columns[current_program] <- final_column
}

marker_availability <- rbindlist(
  marker_availability_list,
  use.names = TRUE,
  fill = TRUE
)

marker_availability_summary <- marker_availability[
  ,
  .(
    requested_markers = .N,
    available_markers = sum(available),
    genes_used = paste(gene[available], collapse = "; ")
  ),
  by = program
]

marker_availability_summary[
  ,
  score_column := unname(score_columns[program])
]

fwrite(
  marker_availability,
  file.path(table3_dir, "PostRT_marker_availability.csv")
)

fwrite(
  marker_availability_summary,
  file.path(table3_dir, "PostRT_marker_availability_summary.csv")
)

required_programs <- c(
  functional_definitions$program,
  cancer_programs
)

missing_programs <- setdiff(
  required_programs,
  names(score_columns)
)

if (length(missing_programs) > 0L) {
  warning(
    "The following programs were not scored and may be omitted: ",
    paste(missing_programs, collapse = ", ")
  )
}

## ------------------------------------------------------------
## 6. Spatial metadata and nearest-cancer distances
## ------------------------------------------------------------

spatial_data <- as.data.table(
  full_object[[]],
  keep.rownames = "cell"
)

spatial_data[, sample_id := as.character(sample_id)]
spatial_data[, fov := as.character(fov)]
spatial_data[, Final_celltype := as.character(Final_celltype)]
spatial_data[, x_px := as.numeric(CenterX_local_px)]
spatial_data[, y_px := as.numeric(CenterY_local_px)]
spatial_data[, x_um := x_px * pixel_size_um]
spatial_data[, y_um := y_px * pixel_size_um]

spatial_data <- spatial_data[
  sample_id == analysis_sample &
    fov %in% target_fovs &
    is.finite(x_px) &
    is.finite(y_px)
]

spatial_data[
  ,
  Region := factor(
    unname(region_mapping[fov]),
    levels = region_levels
  )
]

if (nrow(spatial_data) == 0L) {
  stop("No cells were found for the selected sample and FOVs.")
}

if (any(is.na(spatial_data$Region))) {
  stop("Some cells could not be assigned to a response region.")
}

distance_result_list <- vector(
  mode = "list",
  length = length(target_fovs)
)
names(distance_result_list) <- target_fovs

for (current_fov in target_fovs) {
  current_data <- spatial_data[fov == current_fov]
  cancer_data_fov <- current_data[Final_celltype == "Cancer"]
  noncancer_data_fov <- current_data[Final_celltype != "Cancer"]

  if (nrow(cancer_data_fov) == 0L) {
    stop("No cancer cells were found in FOV ", current_fov, ".")
  }

  if (nrow(noncancer_data_fov) == 0L) {
    stop("No non-cancer cells were found in FOV ", current_fov, ".")
  }

  nearest_cancer <- FNN::get.knnx(
    data = as.matrix(cancer_data_fov[, .(x_px, y_px)]),
    query = as.matrix(noncancer_data_fov[, .(x_px, y_px)]),
    k = 1,
    algorithm = "kd_tree"
  )

  noncancer_data_fov[
    ,
    nearest_cancer_cell := cancer_data_fov$cell[
      as.integer(nearest_cancer$nn.index[, 1])
    ]
  ]

  noncancer_data_fov[
    ,
    distance_to_nearest_cancer_um := as.numeric(
      nearest_cancer$nn.dist[, 1]
    ) * pixel_size_um
  ]

  noncancer_data_fov[
    ,
    tumor_adjacent := distance_to_nearest_cancer_um <= primary_radius_um
  ]

  distance_result_list[[current_fov]] <- noncancer_data_fov
}

noncancer_distance_data <- unique(
  rbindlist(
    distance_result_list,
    use.names = TRUE,
    fill = TRUE
  ),
  by = "cell"
)

fwrite(
  noncancer_distance_data[
    ,
    .(
      cell,
      sample_id,
      fov,
      Region,
      Final_celltype,
      nearest_cancer_cell,
      distance_to_nearest_cancer_um,
      tumor_adjacent
    )
  ],
  file.path(table3_dir, "Figure3_cell_to_nearest_cancer_distance.csv")
)

## ------------------------------------------------------------
## 7. Cell-type-specific tumor-adjacent functional data
## ------------------------------------------------------------

functional_definitions_available <- functional_definitions[
  program %in% names(score_columns)
]

functional_data_list <- vector(
  mode = "list",
  length = nrow(functional_definitions_available)
)

for (definition_index in seq_len(nrow(functional_definitions_available))) {
  current_definition <- functional_definitions_available[definition_index]
  current_score_column <- unname(
    score_columns[current_definition$program]
  )

  current_data <- copy(
    noncancer_distance_data[
      Final_celltype == current_definition$target_celltype
    ]
  )

  if (nrow(current_data) == 0L) {
    next
  }

  current_data[
    ,
    `:=`(
      program = current_definition$program,
      program_display = current_definition$display_label,
      score = as.numeric(get(current_score_column))
    )
  ]

  current_data <- current_data[is.finite(score)]
  current_data[, score_z := standardize_vector(score)]
  functional_data_list[[definition_index]] <- current_data
}

functional_data <- rbindlist(
  functional_data_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(functional_data) == 0L) {
  stop("No cell-type-specific functional-score data were generated.")
}

adjacent_functional_data <- functional_data[tumor_adjacent == TRUE]

functional_statistics <- adjacent_functional_data[
  ,
  {
    nonresponse_values <- score_z[
      as.character(Region) == "Non-response region"
    ]
    response_values <- score_z[
      as.character(Region) == "Response region"
    ]

    .(
      nonresponse_cells = as.integer(length(nonresponse_values)),
      response_cells = as.integer(length(response_values)),
      nonresponse_mean_z = as.numeric(mean(nonresponse_values, na.rm = TRUE)),
      response_mean_z = as.numeric(mean(response_values, na.rm = TRUE)),
      nonresponse_median_z = as.numeric(
        median(nonresponse_values, na.rm = TRUE)
      ),
      response_median_z = as.numeric(
        median(response_values, na.rm = TRUE)
      ),
      response_minus_nonresponse = as.numeric(
        median(response_values, na.rm = TRUE) -
          median(nonresponse_values, na.rm = TRUE)
      ),
      raw_p = safe_wilcoxon(
        nonresponse_values,
        response_values
      )
    )
  },
  by = .(
    program,
    program_display,
    Final_celltype
  )
]

functional_statistics[, FDR := p.adjust(raw_p, method = "BH")]
functional_statistics[, raw_p_label := vapply(raw_p, format_raw_p, character(1))]

fwrite(
  adjacent_functional_data[
    ,
    .(
      cell,
      sample_id,
      fov,
      Region,
      Final_celltype,
      distance_to_nearest_cancer_um,
      program,
      program_display,
      score,
      score_z
    )
  ],
  file.path(table3_dir, "Figure3_tumor_adjacent_functional_scores.csv")
)

fwrite(
  functional_statistics,
  file.path(table3_dir, "Figure3_functional_score_statistics.csv")
)

## ------------------------------------------------------------
## 8. Cancer-cell program data and regional tests
## ------------------------------------------------------------

cancer_data <- copy(
  spatial_data[Final_celltype == "Cancer"]
)

if (nrow(cancer_data) == 0L) {
  stop("No cancer cells were found in the selected FOVs.")
}

available_cancer_programs <- intersect(
  cancer_programs,
  names(score_columns)
)

cancer_score_list <- vector(
  mode = "list",
  length = length(available_cancer_programs)
)

for (program_index in seq_along(available_cancer_programs)) {
  current_program <- available_cancer_programs[program_index]
  current_score_column <- unname(score_columns[current_program])

  current_data <- copy(cancer_data)
  current_data[
    ,
    `:=`(
      program = current_program,
      program_display = unname(
        cancer_program_display[current_program]
      ),
      score = as.numeric(get(current_score_column))
    )
  ]

  current_data <- current_data[is.finite(score)]
  current_data[, score_z := standardize_vector(score)]
  cancer_score_list[[program_index]] <- current_data
}

cancer_score_data <- rbindlist(
  cancer_score_list,
  use.names = TRUE,
  fill = TRUE
)

cancer_score_statistics <- cancer_score_data[
  ,
  {
    nonresponse_values <- score_z[
      as.character(Region) == "Non-response region"
    ]
    response_values <- score_z[
      as.character(Region) == "Response region"
    ]

    .(
      nonresponse_cells = as.integer(length(nonresponse_values)),
      response_cells = as.integer(length(response_values)),
      nonresponse_median_z = as.numeric(
        median(nonresponse_values, na.rm = TRUE)
      ),
      response_median_z = as.numeric(
        median(response_values, na.rm = TRUE)
      ),
      response_minus_nonresponse = as.numeric(
        median(response_values, na.rm = TRUE) -
          median(nonresponse_values, na.rm = TRUE)
      ),
      raw_p = safe_wilcoxon(
        nonresponse_values,
        response_values
      )
    )
  },
  by = .(
    program,
    program_display
  )
]

cancer_score_statistics[, FDR := p.adjust(raw_p, method = "BH")]
cancer_score_statistics[, raw_p_label := vapply(raw_p, format_raw_p, character(1))]

fwrite(
  cancer_score_data[
    ,
    .(
      cell,
      sample_id,
      fov,
      Region,
      program,
      program_display,
      score,
      score_z
    )
  ],
  file.path(table4_dir, "Figure4_cancer_program_cell_scores.csv")
)

fwrite(
  cancer_score_statistics,
  file.path(table4_dir, "Figure4_cancer_program_statistics.csv")
)

## ------------------------------------------------------------
## 9. Selected cancer-cell gene effect sizes
## ------------------------------------------------------------

normalized_layer <- LayerData(
  object = full_object,
  assay = "Nanostring",
  layer = "data"
)

selected_gene_table <- rbindlist(
  lapply(
    names(selected_gene_groups),
    function(current_group) {
      data.table(
        gene_group = current_group,
        gene = selected_gene_groups[[current_group]]
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

selected_gene_table <- selected_gene_table[
  gene %in% rownames(normalized_layer)
]

selected_gene_table <- unique(selected_gene_table, by = "gene")

if (nrow(selected_gene_table) == 0L) {
  stop("None of the selected Figure 4 genes are in the data layer.")
}

cancer_cells <- intersect(
  cancer_data$cell,
  colnames(normalized_layer)
)

gene_expression_matrix <- as.matrix(
  normalized_layer[
    selected_gene_table$gene,
    cancer_cells,
    drop = FALSE
  ]
)

cancer_region_lookup <- setNames(
  as.character(cancer_data$Region),
  cancer_data$cell
)

nonresponse_cells <- cancer_cells[
  cancer_region_lookup[cancer_cells] == "Non-response region"
]

response_cells <- cancer_cells[
  cancer_region_lookup[cancer_cells] == "Response region"
]

if (length(nonresponse_cells) < 2L || length(response_cells) < 2L) {
  stop("Insufficient cancer cells for the regional gene comparison.")
}

gene_effect_list <- vector(
  mode = "list",
  length = nrow(gene_expression_matrix)
)

for (gene_index in seq_len(nrow(gene_expression_matrix))) {
  current_gene <- rownames(gene_expression_matrix)[gene_index]
  nonresponse_values <- as.numeric(
    gene_expression_matrix[gene_index, nonresponse_cells]
  )
  response_values <- as.numeric(
    gene_expression_matrix[gene_index, response_cells]
  )

  mean_nonresponse_linear <- mean(
    expm1(nonresponse_values),
    na.rm = TRUE
  )

  mean_response_linear <- mean(
    expm1(response_values),
    na.rm = TRUE
  )

  gene_effect_list[[gene_index]] <- data.table(
    gene = current_gene,
    mean_nonresponse = as.numeric(mean_nonresponse_linear),
    mean_response = as.numeric(mean_response_linear),
    pct_nonresponse = as.numeric(mean(nonresponse_values > 0) * 100),
    pct_response = as.numeric(mean(response_values > 0) * 100),
    log2FC_response_vs_nonresponse = as.numeric(
      log2(
        (mean_response_linear + 1e-6) /
          (mean_nonresponse_linear + 1e-6)
      )
    ),
    raw_p = safe_wilcoxon(
      nonresponse_values,
      response_values
    )
  )
}

selected_gene_effects <- rbindlist(
  gene_effect_list,
  use.names = TRUE,
  fill = TRUE
)

selected_gene_effects <- merge(
  selected_gene_effects,
  selected_gene_table,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

selected_gene_effects[, FDR := p.adjust(raw_p, method = "BH")]
selected_gene_effects[, significance := vapply(raw_p, raw_p_symbol, character(1))]
selected_gene_effects[
  ,
  enriched_region := factor(
    fifelse(
      log2FC_response_vs_nonresponse >= 0,
      "Response region",
      "Non-response region"
    ),
    levels = region_levels
  )
]

selected_gene_effects[
  ,
  gene_group := factor(
    gene_group,
    levels = names(selected_gene_groups)
  )
]

selected_gene_effects[
  ,
  gene_order := match(
    gene,
    unlist(selected_gene_groups, use.names = FALSE)
  )
]

selected_gene_effects <- selected_gene_effects[order(gene_order)]
selected_gene_effects[
  ,
  gene_plot := factor(
    gene,
    levels = rev(unique(gene))
  )
]

fwrite(
  selected_gene_effects,
  file.path(table4_dir, "Figure4_selected_cancer_gene_effects.csv")
)

## ------------------------------------------------------------
## 10. Cancer IFN signaling and local ECM-CAF burden for Figure 4F
## ------------------------------------------------------------

ecm_caf_spatial_data <- copy(
  functional_data[
    program == "CAF_ECM_remodeling"
  ]
)

cancer_ifn_spatial_data <- copy(
  cancer_score_data[
    program == "Cancer_IFN_adaptive_resistance"
  ]
)

if (nrow(ecm_caf_spatial_data) == 0L) {
  stop("No ECM-CAF cells are available for the local coupling analysis.")
}

if (nrow(cancer_ifn_spatial_data) == 0L) {
  stop("No cancer IFN scores are available for the local coupling analysis.")
}

setnames(
  cancer_ifn_spatial_data,
  old = c("score", "score_z"),
  new = c("cancer_ifn_score", "cancer_ifn_score_z")
)

ecm_ifn_coupling_list <- vector(
  mode = "list",
  length = length(region_levels)
)

names(ecm_ifn_coupling_list) <- region_levels

for (current_region in region_levels) {
  current_cancer <- copy(
    cancer_ifn_spatial_data[
      as.character(Region) == current_region
    ]
  )

  current_ecm_caf <- copy(
    ecm_caf_spatial_data[
      as.character(Region) == current_region
    ]
  )

  current_cancer[, local_ecm_caf_count := 0L]
  current_cancer[, local_mean_ecm_score_z := NA_real_]
  current_cancer[, local_median_ecm_score_z := NA_real_]

  if (nrow(current_cancer) == 0L) {
    ecm_ifn_coupling_list[[current_region]] <- current_cancer
    next
  }

  if (nrow(current_ecm_caf) > 0L) {
    neighbor_result <- FNN::get.knnx(
      data = as.matrix(current_ecm_caf[, .(x_um, y_um)]),
      query = as.matrix(current_cancer[, .(x_um, y_um)]),
      k = nrow(current_ecm_caf),
      algorithm = "kd_tree"
    )

    neighbor_indices <- neighbor_result$nn.index
    neighbor_distances <- neighbor_result$nn.dist

    if (is.null(dim(neighbor_indices))) {
      neighbor_indices <- matrix(neighbor_indices, ncol = 1L)
      neighbor_distances <- matrix(neighbor_distances, ncol = 1L)
    }

    local_counts <- integer(nrow(current_cancer))
    local_mean_scores <- rep(NA_real_, nrow(current_cancer))
    local_median_scores <- rep(NA_real_, nrow(current_cancer))

    for (cancer_index in seq_len(nrow(current_cancer))) {
      inside_radius <- which(
        neighbor_distances[cancer_index, ] <= primary_radius_um
      )

      local_counts[cancer_index] <- length(inside_radius)

      if (length(inside_radius) > 0L) {
        ecm_rows <- neighbor_indices[
          cancer_index,
          inside_radius,
          drop = TRUE
        ]

        local_scores <- as.numeric(
          current_ecm_caf$score_z[ecm_rows]
        )

        local_mean_scores[cancer_index] <- mean(
          local_scores,
          na.rm = TRUE
        )

        local_median_scores[cancer_index] <- median(
          local_scores,
          na.rm = TRUE
        )
      }
    }

    current_cancer[
      ,
      local_ecm_caf_count := as.integer(local_counts)
    ]

    current_cancer[
      ,
      local_mean_ecm_score_z := as.numeric(local_mean_scores)
    ]

    current_cancer[
      ,
      local_median_ecm_score_z := as.numeric(local_median_scores)
    ]
  }

  current_cancer[
    ,
    log1p_local_ecm_caf_count := log1p(local_ecm_caf_count)
  ]

  ecm_ifn_coupling_list[[current_region]] <- current_cancer
}

ecm_ifn_coupling_data <- rbindlist(
  ecm_ifn_coupling_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(ecm_ifn_coupling_data) == 0L) {
  stop("No cancer cells were available for the ECM-CAF coupling analysis.")
}

ecm_ifn_coupling_data[
  ,
  Region := factor(as.character(Region), levels = region_levels)
]

ecm_ifn_correlation_statistics <- ecm_ifn_coupling_data[
  ,
  {
    count_test <- safe_spearman(
      local_ecm_caf_count,
      cancer_ifn_score_z
    )

    activity_rows <-
      local_ecm_caf_count > 0L &
      is.finite(local_mean_ecm_score_z)

    activity_test <- safe_spearman(
      local_mean_ecm_score_z[activity_rows],
      cancer_ifn_score_z[activity_rows]
    )

    .(
      cancer_cells = as.integer(.N),
      cancer_cells_with_local_ecm_caf = as.integer(sum(activity_rows)),
      fraction_with_local_ecm_caf = as.numeric(mean(activity_rows)),
      median_local_ecm_caf_count = as.numeric(
        median(local_ecm_caf_count, na.rm = TRUE)
      ),
      rho_ifn_vs_ecm_caf_count = as.numeric(count_test$rho),
      raw_p_ifn_vs_ecm_caf_count = as.numeric(count_test$raw_p),
      rho_ifn_vs_mean_ecm_activity = as.numeric(activity_test$rho),
      raw_p_ifn_vs_mean_ecm_activity = as.numeric(activity_test$raw_p)
    )
  },
  by = Region
]

ecm_ifn_correlation_statistics[
  ,
  count_correlation_label := paste0(
    "Spearman ρ = ",
    ifelse(
      is.finite(rho_ifn_vs_ecm_caf_count),
      formatC(
        rho_ifn_vs_ecm_caf_count,
        format = "f",
        digits = 2
      ),
      "NA"
    ),
    "\n",
    vapply(
      raw_p_ifn_vs_ecm_caf_count,
      format_raw_p,
      character(1)
    )
  )
]

ecm_ifn_adjusted_model <- lm(
  cancer_ifn_score_z ~
    log1p_local_ecm_caf_count * Region,
  data = ecm_ifn_coupling_data
)

ecm_ifn_adjusted_model_table <- as.data.table(
  summary(ecm_ifn_adjusted_model)$coefficients,
  keep.rownames = "term"
)

setnames(
  ecm_ifn_adjusted_model_table,
  old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
  new = c("estimate", "standard_error", "t_value", "raw_p")
)

fwrite(
  ecm_ifn_coupling_data[
    ,
    .(
      cell,
      sample_id,
      fov,
      Region,
      x_um,
      y_um,
      cancer_ifn_score,
      cancer_ifn_score_z,
      local_ecm_caf_count,
      local_mean_ecm_score_z,
      local_median_ecm_score_z,
      log1p_local_ecm_caf_count
    )
  ],
  file.path(
    table4_dir,
    "Figure4_cancer_IFN_local_ECM_CAF_data.csv"
  )
)

fwrite(
  ecm_ifn_correlation_statistics,
  file.path(
    table4_dir,
    "Figure4_cancer_IFN_local_ECM_CAF_correlations.csv"
  )
)

fwrite(
  ecm_ifn_adjusted_model_table,
  file.path(
    table4_dir,
    "Figure4_cancer_IFN_local_ECM_CAF_adjusted_model.csv"
  )
)

## ------------------------------------------------------------
## 11. Selected functional genes for Figure 3I-J and Supplementary Figure 3E-F
## ------------------------------------------------------------

figure3_gene_sets <- list(
  "ECM-CAF remodeling" = c(
    "COL1A1", "COL1A2", "COL3A1", "FN1",
    "THBS2", "MMP1", "VCAN", "INHBA"
  ),
  "Myeloid inflammation" = c(
    "IL1B", "CXCL8", "S100A8", "S100A9",
    "IL1RN", "CCL3/L1/L3", "PTGS2", "CSF3R"
  ),
  "T-cell effector/activation" = c(
    "CD69", "ICOS", "TNFRSF9", "IFNG",
    "GZMK", "CCL5", "IL32", "RUNX3"
  ),
  "Cytotoxic NK/T activity" = c(
    "NKG7", "GNLY", "PRF1", "GZMB", "GZMH",
    "CCL5", "KLRD1", "KLRK1", "EOMES"
  )
)

figure3_program_mapping <- c(
  "ECM-CAF remodeling" = "CAF_ECM_remodeling",
  "Myeloid inflammation" = "Myeloid_inflammatory",
  "T-cell effector/activation" = "Tcell_effector_activation",
  "Cytotoxic NK/T activity" = "Cytotoxic_NKT_cytotoxicity"
)

figure3_gene_table <- rbindlist(
  lapply(
    names(figure3_gene_sets),
    function(current_group) {
      data.table(
        gene_group = current_group,
        program = unname(
          figure3_program_mapping[current_group]
        ),
        gene = figure3_gene_sets[[current_group]]
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

figure3_gene_table <- figure3_gene_table[
  gene %in% rownames(normalized_layer)
]

if (nrow(figure3_gene_table) == 0L) {
  stop("None of the selected Figure 3 genes are in the data layer.")
}

figure3_gene_effect_list <- vector(
  mode = "list",
  length = nrow(figure3_gene_table)
)

for (gene_index in seq_len(nrow(figure3_gene_table))) {
  current_definition <- figure3_gene_table[gene_index]

  current_cell_data <- adjacent_functional_data[
    program == current_definition$program
  ]

  current_cells <- intersect(
    current_cell_data$cell,
    colnames(normalized_layer)
  )

  current_region_lookup <- setNames(
    as.character(current_cell_data$Region),
    current_cell_data$cell
  )

  current_expression <- as.numeric(
    normalized_layer[
      current_definition$gene,
      current_cells,
      drop = TRUE
    ]
  )

  nonresponse_values <- current_expression[
    current_region_lookup[current_cells] ==
      "Non-response region"
  ]

  response_values <- current_expression[
    current_region_lookup[current_cells] ==
      "Response region"
  ]

  mean_nonresponse_linear <- mean(
    expm1(nonresponse_values),
    na.rm = TRUE
  )

  mean_response_linear <- mean(
    expm1(response_values),
    na.rm = TRUE
  )

  figure3_gene_effect_list[[gene_index]] <- data.table(
    gene_group = current_definition$gene_group,
    program = current_definition$program,
    gene = current_definition$gene,
    mean_nonresponse = as.numeric(mean_nonresponse_linear),
    mean_response = as.numeric(mean_response_linear),
    pct_nonresponse = as.numeric(
      mean(nonresponse_values > 0) * 100
    ),
    pct_response = as.numeric(
      mean(response_values > 0) * 100
    ),
    log2FC_response_vs_nonresponse = as.numeric(
      log2(
        (mean_response_linear + 1e-6) /
          (mean_nonresponse_linear + 1e-6)
      )
    ),
    raw_p = safe_wilcoxon(
      nonresponse_values,
      response_values
    )
  )
}

figure3_selected_gene_effects <- rbindlist(
  figure3_gene_effect_list,
  use.names = TRUE,
  fill = TRUE
)

figure3_selected_gene_effects[
  ,
  FDR := p.adjust(raw_p, method = "BH")
]
figure3_selected_gene_effects[
  ,
  significance := vapply(raw_p, raw_p_symbol, character(1))
]
figure3_selected_gene_effects[
  ,
  enriched_region := factor(
    fifelse(
      log2FC_response_vs_nonresponse >= 0,
      "Response region",
      "Non-response region"
    ),
    levels = region_levels
  )
]
figure3_selected_gene_effects[
  ,
  gene_group := factor(
    gene_group,
    levels = names(figure3_gene_sets)
  )
]
figure3_selected_gene_effects[
  ,
  gene_order := match(
    gene,
    unlist(figure3_gene_sets, use.names = FALSE)
  )
]

figure3_selected_gene_effects <-
  figure3_selected_gene_effects[order(gene_order)]

figure3_selected_gene_effects[
  ,
  gene_plot := factor(
    gene,
    levels = rev(unique(gene))
  )
]

fwrite(
  figure3_selected_gene_effects,
  file.path(table3_dir, "Figure3_selected_functional_gene_effects.csv")
)

## ------------------------------------------------------------
## 12. Cancer-centered spatial ligand-receptor analysis
## ------------------------------------------------------------

lymphoid_celltypes_text <- paste(
  c(
    "Imm_Tcell",
    "Imm_Cytotoxic_NK/T",
    "Imm_NK"
  ),
  collapse = ";"
)

interaction_candidates <- data.table(
  interaction_class = c(
    rep("Cancer to Myeloid", 2),
    rep("Cancer to T/NK", 7),
    rep("ECM-CAF to Cancer", 3),
    rep("Myeloid to Cancer", 3),
    "T/NK to Cancer"
  ),
  sender_group = c(
    rep("Cancer", 9),
    rep("CAF_ECM", 3),
    rep("Imm_Myeloid", 3),
    "T/NK"
  ),
  sender_celltypes = c(
    rep("Cancer", 9),
    rep("CAF_ECM", 3),
    rep("Imm_Myeloid", 3),
    lymphoid_celltypes_text
  ),
  ligand = c(
    "CSF1", "CD47",
    "CXCL9", "CXCL10", "CD274", "PVR", "NECTIN2",
    "LGALS9", "ICAM1",
    "FN1", "THBS1", "TGFB1",
    "IL1B", "TGFB1", "SPP1",
    "IFNG"
  ),
  receptor = c(
    "CSF1R", "SIRPA",
    "CXCR3", "CXCR3", "PDCD1", "TIGIT", "TIGIT",
    "HAVCR2", "ITGAL",
    "CD44", "CD47", "TGFBR2",
    "IL1R1", "TGFBR2", "CD44",
    "IFNGR1"
  ),
  receiver_group = c(
    rep("Imm_Myeloid", 2),
    rep("T/NK", 7),
    rep("Cancer", 7)
  ),
  receiver_celltypes = c(
    rep("Imm_Myeloid", 2),
    rep(lymphoid_celltypes_text, 7),
    rep("Cancer", 7)
  )
)

interaction_candidates[, priority := seq_len(.N)]
interaction_candidates[
  ,
  interaction_id := paste(
    sender_group,
    ligand,
    receptor,
    receiver_group,
    sep = "__"
  )
]
interaction_candidates[
  ,
  pair_label := paste0(ligand, " → ", receptor)
]
interaction_candidates[
  ,
  ligand_available := ligand %in% rownames(normalized_layer)
]
interaction_candidates[
  ,
  receptor_available := receptor %in% rownames(normalized_layer)
]
interaction_candidates[
  ,
  genes_available := ligand_available & receptor_available
]

available_interactions <- interaction_candidates[
  genes_available == TRUE
]

if (nrow(available_interactions) == 0L) {
  stop(
    "None of the prespecified ligand-receptor pairs are available ",
    "in the Nanostring panel."
  )
}

interaction_score_list <- vector(
  mode = "list",
  length = nrow(available_interactions) * length(region_levels)
)

interaction_score_index <- 0L

for (interaction_index in seq_len(nrow(available_interactions))) {
  current_interaction <- available_interactions[interaction_index]

  current_sender_types <- strsplit(
    current_interaction$sender_celltypes,
    split = ";",
    fixed = TRUE
  )[[1]]

  current_receiver_types <- strsplit(
    current_interaction$receiver_celltypes,
    split = ";",
    fixed = TRUE
  )[[1]]

  for (current_region in region_levels) {
    interaction_score_index <- interaction_score_index + 1L

    current_sender <- spatial_data[
      as.character(Region) == current_region &
        Final_celltype %in% current_sender_types
    ]

    current_receiver <- spatial_data[
      as.character(Region) == current_region &
        Final_celltype %in% current_receiver_types
    ]

    if (nrow(current_sender) == 0L) {
      next
    }

    current_scores <- calculate_sender_local_interaction(
      sender_table = current_sender,
      receiver_table = current_receiver,
      expression_matrix = normalized_layer,
      ligand_gene = current_interaction$ligand,
      receptor_gene = current_interaction$receptor,
      radius_um = primary_radius_um
    )

    current_scores[
      ,
      `:=`(
        interaction_id = current_interaction$interaction_id,
        interaction_class = current_interaction$interaction_class,
        sender_group = current_interaction$sender_group,
        receiver_group = current_interaction$receiver_group,
        ligand = current_interaction$ligand,
        receptor = current_interaction$receptor,
        pair_label = current_interaction$pair_label,
        priority = current_interaction$priority,
        sender_cells_in_region = as.integer(nrow(current_sender)),
        receiver_cells_in_region = as.integer(nrow(current_receiver))
      )
    ]

    interaction_score_list[[interaction_score_index]] <- current_scores
  }
}

interaction_sender_scores <- rbindlist(
  interaction_score_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(interaction_sender_scores) == 0L) {
  stop("No spatial ligand-receptor scores could be calculated.")
}

interaction_statistics <- interaction_sender_scores[
  ,
  {
    nonresponse_values <- local_interaction_score[
      as.character(Region) == "Non-response region"
    ]

    response_values <- local_interaction_score[
      as.character(Region) == "Response region"
    ]

    nonresponse_codetection <- ligand_receptor_codetected[
      as.character(Region) == "Non-response region"
    ]

    response_codetection <- ligand_receptor_codetected[
      as.character(Region) == "Response region"
    ]

    nonresponse_sender_count <- unique(
      sender_cells_in_region[
        as.character(Region) == "Non-response region"
      ]
    )

    response_sender_count <- unique(
      sender_cells_in_region[
        as.character(Region) == "Response region"
      ]
    )

    nonresponse_receiver_count <- unique(
      receiver_cells_in_region[
        as.character(Region) == "Non-response region"
      ]
    )

    response_receiver_count <- unique(
      receiver_cells_in_region[
        as.character(Region) == "Response region"
      ]
    )

    mean_nonresponse <- mean(
      nonresponse_values,
      na.rm = TRUE
    )

    mean_response <- mean(
      response_values,
      na.rm = TRUE
    )

    .(
      nonresponse_sender_cells = as.integer(
        ifelse(length(nonresponse_sender_count) == 0L, 0L,
               nonresponse_sender_count[1])
      ),
      response_sender_cells = as.integer(
        ifelse(length(response_sender_count) == 0L, 0L,
               response_sender_count[1])
      ),
      nonresponse_receiver_cells = as.integer(
        ifelse(length(nonresponse_receiver_count) == 0L, 0L,
               nonresponse_receiver_count[1])
      ),
      response_receiver_cells = as.integer(
        ifelse(length(response_receiver_count) == 0L, 0L,
               response_receiver_count[1])
      ),
      mean_score_nonresponse = as.numeric(mean_nonresponse),
      mean_score_response = as.numeric(mean_response),
      median_score_nonresponse = as.numeric(
        median(nonresponse_values, na.rm = TRUE)
      ),
      median_score_response = as.numeric(
        median(response_values, na.rm = TRUE)
      ),
      codetection_nonresponse = as.numeric(
        mean(nonresponse_codetection, na.rm = TRUE)
      ),
      codetection_response = as.numeric(
        mean(response_codetection, na.rm = TRUE)
      ),
      log2FC_response_vs_nonresponse = as.numeric(
        log2(
          (mean_response + interaction_pseudocount) /
            (mean_nonresponse + interaction_pseudocount)
        )
      ),
      rank_biserial_response_vs_nonresponse =
        safe_rank_biserial(
          response_values = response_values,
          nonresponse_values = nonresponse_values
        ),
      raw_p = safe_wilcoxon(
        nonresponse_values,
        response_values
      )
    )
  },
  by = .(
    interaction_id,
    interaction_class,
    sender_group,
    receiver_group,
    ligand,
    receptor,
    pair_label,
    priority
  )
]

interaction_statistics[, FDR := p.adjust(raw_p, method = "BH")]
interaction_statistics[
  ,
  significance := vapply(raw_p, raw_p_symbol, character(1))
]
interaction_statistics[
  ,
  maximum_codetection := pmax(
    codetection_nonresponse,
    codetection_response,
    na.rm = TRUE
  )
]
interaction_statistics[
  ,
  eligible_for_main_figure :=
    nonresponse_sender_cells >= minimum_interaction_sender_cells &
    response_sender_cells >= minimum_interaction_sender_cells &
    nonresponse_receiver_cells >= minimum_interaction_receiver_cells &
    response_receiver_cells >= minimum_interaction_receiver_cells &
    maximum_codetection >= minimum_codetection_fraction &
    is.finite(rank_biserial_response_vs_nonresponse)
]

interaction_plot_data <- interaction_statistics[
  eligible_for_main_figure == TRUE
]

if (nrow(interaction_plot_data) == 0L) {
  warning(
    "No interaction passed the co-detection threshold; ",
    "the threshold was relaxed for plotting."
  )

  interaction_plot_data <- interaction_statistics[
    nonresponse_sender_cells >= minimum_interaction_sender_cells &
      response_sender_cells >= minimum_interaction_sender_cells &
      nonresponse_receiver_cells >= minimum_interaction_receiver_cells &
      response_receiver_cells >= minimum_interaction_receiver_cells &
      is.finite(rank_biserial_response_vs_nonresponse)
  ]
}

if (nrow(interaction_plot_data) == 0L) {
  stop(
    "No ligand-receptor interaction was eligible for ",
    "Supplementary Figure 4."
  )
}

interaction_plot_data[
  ,
  absolute_effect := abs(
    rank_biserial_response_vs_nonresponse
  )
]

setorder(
  interaction_plot_data,
  -absolute_effect,
  priority
)

interaction_plot_data <- head(
  interaction_plot_data,
  maximum_main_interactions
)

interaction_plot_data[
  ,
  enriched_region := factor(
    fifelse(
      rank_biserial_response_vs_nonresponse >= 0,
      "Response region",
      "Non-response region"
    ),
    levels = region_levels
  )
]
interaction_plot_data[
  ,
  interaction_class := factor(
    interaction_class,
    levels = unique(interaction_candidates$interaction_class)
  )
]
interaction_plot_data[
  ,
  pair_plot := factor(
    pair_label,
    levels = rev(unique(pair_label))
  )
]
interaction_plot_data[
  ,
  codetection_percent := 100 * maximum_codetection
]
interaction_plot_data[
  ,
  significance_x := rank_biserial_response_vs_nonresponse +
    fifelse(
      rank_biserial_response_vs_nonresponse >= 0,
      0.045,
      -0.045
    )
]

fwrite(
  interaction_candidates,
  file.path(
    table4_dir,
    "SupplementaryFigure4_ligand_receptor_panel_coverage.csv"
  )
)

fwrite(
  interaction_sender_scores[
    ,
    .(
      cell,
      sample_id,
      fov,
      Region,
      Final_celltype,
      interaction_id,
      interaction_class,
      sender_group,
      receiver_group,
      ligand,
      receptor,
      ligand_expression,
      local_receiver_count,
      local_mean_receptor_expression,
      local_receptor_positive_fraction,
      local_interaction_score,
      ligand_receptor_codetected
    )
  ],
  file.path(
    table4_dir,
    "SupplementaryFigure4_sender_local_interaction_scores.csv"
  )
)

fwrite(
  interaction_statistics,
  file.path(
    table4_dir,
    "SupplementaryFigure4_spatial_interaction_statistics.csv"
  )
)

## ------------------------------------------------------------
## 13. Plot helper functions
## ------------------------------------------------------------

make_combined_tme_spatial_plot <- function(
  region_name,
  plot_title,
  show_legend = TRUE,
  maximum_distance_um = primary_radius_um
) {
  current_background <- spatial_data[
    as.character(Region) == region_name
  ]

  current_cancer <- current_background[
    Final_celltype == "Cancer"
  ]

  current_target <- copy(
    functional_data[
      program %in% c(
        "CAF_ECM_remodeling",
        "Myeloid_inflammatory"
      ) &
        as.character(Region) == region_name &
        distance_to_nearest_cancer_um <= maximum_distance_um
    ]
  )

  current_target[
    ,
    functional_population := factor(
      fifelse(
        program == "CAF_ECM_remodeling",
        "ECM-CAF",
        "Inflammatory myeloid"
      ),
      levels = c(
        "Cancer",
        "ECM-CAF",
        "Inflammatory myeloid"
      )
    )
  ]

  target_cells <- unique(current_target$cell)

  current_background_grey <- current_background[
    !cell %in% target_cells &
      Final_celltype != "Cancer"
  ]

  common_limit <- score_limit(
    functional_data[
      program %in% c(
        "CAF_ECM_remodeling",
        "Myeloid_inflammatory"
      ) &
        distance_to_nearest_cancer_um <= maximum_distance_um,
      score_z
    ]
  )

  current_plot <- ggplot() +
    geom_point(
      data = current_background_grey,
      aes(x = x_um, y = y_um),
      colour = "grey88",
      size = 0.32,
      alpha = 0.40
    ) +
    geom_point(
      data = current_cancer,
      aes(
        x = x_um,
        y = y_um,
        shape = "Cancer"
      ),
      colour = unname(celltype_colors[["Cancer"]]),
      fill = unname(celltype_colors[["Cancer"]]),
      size = 0.58,
      alpha = 0.88
    ) +
    geom_point(
      data = current_target,
      aes(
        x = x_um,
        y = y_um,
        fill = score_z,
        shape = functional_population,
        colour = functional_population
      ),
      stroke = 0.18,
      size = 1.52,
      alpha = 0.96
    ) +
    scale_fill_gradient2(
      low = score_low_color,
      mid = "white",
      high = score_high_color,
      midpoint = 0,
      limits = c(-common_limit, common_limit),
      oob = squish,
      name = "Standardized\nmodule score"
    ) +
    scale_colour_manual(
      values = c(
        "ECM-CAF" = unname(celltype_colors[["CAF_ECM"]]),
        "Inflammatory myeloid" =
          unname(celltype_colors[["Imm_Myeloid"]])
      ),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        "Cancer" = 22,
        "ECM-CAF" = 21,
        "Inflammatory myeloid" = 24
      ),
      breaks = c(
        "Cancer",
        "ECM-CAF",
        "Inflammatory myeloid"
      ),
      name = NULL
    ) +
    coord_equal(clip = "off") +
    labs(
      title = plot_title,
      x = "X coordinate (µm)",
      y = "Y coordinate (µm)",
      caption = paste0(
        "Colored cells: ECM-CAF and myeloid cells within ",
        formatC(
          maximum_distance_um,
          format = "f",
          digits = 2
        ),
        " µm of the nearest cancer cell"
      )
    ) +
    theme_manuscript +
    theme(
      panel.grid = element_blank(),
      legend.text = element_text(
        size = 10.5,
        colour = "black"
      ),
      legend.title = element_text(
        size = 10.5,
        face = "bold",
        colour = "black"
      ),
      legend.key.height = grid::unit(0.48, "cm"),
      legend.key.width = grid::unit(0.48, "cm"),
      plot.caption = element_text(
        hjust = 1,
        size = 7.2,
        colour = "grey30"
      )
    ) +
    guides(
      fill = guide_colourbar(
        barheight = grid::unit(2.5, "cm"),
        barwidth = grid::unit(0.48, "cm")
      ),
      shape = guide_legend(
        override.aes = list(
          size = 3.2,
          alpha = 1,
          fill = c(
            unname(celltype_colors[["Cancer"]]),
            "white",
            "white"
          ),
          colour = c(
            unname(celltype_colors[["Cancer"]]),
            unname(celltype_colors[["CAF_ECM"]]),
            unname(celltype_colors[["Imm_Myeloid"]])
          )
        )
      )
    )

  if (!isTRUE(show_legend)) {
    current_plot <- current_plot + theme(legend.position = "none")
  }

  current_plot
}

make_combined_immune_spatial_plot <- function(
  region_name,
  plot_title,
  show_legend = TRUE
) {
  current_background <- spatial_data[
    as.character(Region) == region_name
  ]

  current_cancer <- current_background[
    Final_celltype == "Cancer"
  ]

  current_target <- copy(
    adjacent_functional_data[
      program %in% c(
        "Tcell_effector_activation",
        "Cytotoxic_NKT_cytotoxicity"
      ) &
        as.character(Region) == region_name
    ]
  )

  current_target[
    ,
    functional_population := factor(
      fifelse(
        program == "Tcell_effector_activation",
        "T cell",
        "Cytotoxic NK/T"
      ),
      levels = c(
        "T cell",
        "Cytotoxic NK/T"
      )
    )
  ]

  target_cells <- unique(current_target$cell)

  current_background_grey <- current_background[
    !cell %in% target_cells &
      Final_celltype != "Cancer"
  ]

  common_limit <- score_limit(
    adjacent_functional_data[
      program %in% c(
        "Tcell_effector_activation",
        "Cytotoxic_NKT_cytotoxicity"
      ),
      score_z
    ]
  )

  current_plot <- ggplot() +
    geom_point(
      data = current_background_grey,
      aes(x = x_um, y = y_um),
      colour = "grey88",
      size = 0.32,
      alpha = 0.40
    ) +
    geom_point(
      data = current_cancer,
      aes(x = x_um, y = y_um),
      colour = unname(celltype_colors[["Cancer"]]),
      size = 0.58,
      alpha = 0.88
    ) +
    geom_point(
      data = current_target,
      aes(
        x = x_um,
        y = y_um,
        fill = score_z,
        shape = functional_population,
        colour = functional_population
      ),
      stroke = 0.18,
      size = 1.65,
      alpha = 0.98
    ) +
    scale_fill_gradient2(
      low = score_low_color,
      mid = "white",
      high = score_high_color,
      midpoint = 0,
      limits = c(-common_limit, common_limit),
      oob = squish,
      name = "Standardized\nmodule score"
    ) +
    scale_colour_manual(
      values = c(
        "T cell" = unname(celltype_colors[["Imm_Tcell"]]),
        "Cytotoxic NK/T" =
          unname(celltype_colors[["Imm_Cytotoxic_NK/T"]])
      ),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        "T cell" = 21,
        "Cytotoxic NK/T" = 24
      ),
      name = NULL
    ) +
    coord_equal(clip = "off") +
    labs(
      title = plot_title,
      x = "X coordinate (µm)",
      y = "Y coordinate (µm)"
    ) +
    theme_manuscript +
    theme(
      panel.grid = element_blank(),
      legend.text = element_text(
        size = 10.5,
        colour = "black"
      ),
      legend.title = element_text(
        size = 10.5,
        face = "bold",
        colour = "black"
      ),
      legend.key.height = grid::unit(0.48, "cm"),
      legend.key.width = grid::unit(0.48, "cm")
    ) +
    guides(
      fill = guide_colourbar(
        barheight = grid::unit(2.5, "cm"),
        barwidth = grid::unit(0.48, "cm")
      ),
      shape = guide_legend(
        override.aes = list(
          size = 3.2,
          alpha = 1,
          colour = c(
            unname(celltype_colors[["Imm_Tcell"]]),
            unname(celltype_colors[["Imm_Cytotoxic_NK/T"]])
          )
        )
      )
    )

  if (!isTRUE(show_legend)) {
    current_plot <- current_plot + theme(legend.position = "none")
  }

  current_plot
}

make_spatial_function_plot <- function(
  program_name,
  region_name,
  plot_title,
  show_legend = TRUE
) {
  current_background <- spatial_data[
    as.character(Region) == region_name
  ]

  current_noncancer <- current_background[
    Final_celltype != "Cancer"
  ]

  current_cancer <- current_background[
    Final_celltype == "Cancer"
  ]

  current_target <- adjacent_functional_data[
    program == program_name &
      as.character(Region) == region_name
  ]

  common_limit <- score_limit(
    adjacent_functional_data[
      program == program_name,
      score_z
    ]
  )

  current_plot <- ggplot() +
    geom_point(
      data = current_noncancer,
      aes(x = x_um, y = y_um),
      colour = "grey86",
      size = 0.36,
      alpha = 0.45
    ) +
    geom_point(
      data = current_cancer,
      aes(x = x_um, y = y_um),
      colour = "black",
      size = 0.62,
      alpha = 0.90
    ) +
    geom_point(
      data = current_target,
      aes(x = x_um, y = y_um, colour = score_z),
      size = 1.45,
      alpha = 0.98
    ) +
    scale_colour_gradient2(
      low = score_low_color,
      mid = "white",
      high = score_high_color,
      midpoint = 0,
      limits = c(-common_limit, common_limit),
      oob = squish,
      name = "Standardized\nmodule score"
    ) +
    coord_equal(clip = "off") +
    labs(
      title = plot_title,
      x = "X coordinate (µm)",
      y = "Y coordinate (µm)"
    ) +
    theme_manuscript +
    theme(panel.grid = element_blank())

  if (!isTRUE(show_legend)) {
    current_plot <- current_plot + theme(legend.position = "none")
  }

  current_plot
}

make_function_violin <- function(
  program_name,
  plot_title
) {
  current_data <- adjacent_functional_data[
    program == program_name
  ]

  current_statistic <- functional_statistics[
    program == program_name
  ]

  if (nrow(current_data) == 0L || nrow(current_statistic) == 0L) {
    stop("No data are available for ", program_name, ".")
  }

  y_maximum <- max(current_data$score_z, na.rm = TRUE)
  y_range <- diff(range(current_data$score_z, na.rm = TRUE))

  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 1
  }

  ggplot(
    current_data,
    aes(x = Region, y = score_z, fill = Region)
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      alpha = 0.55,
      linewidth = 0.4
    ) +
    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.45
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,
      size = 2.8,
      fill = "white",
      colour = "black"
    ) +
    annotate(
      geom = "text",
      x = 1.5,
      y = y_maximum + 0.12 * y_range,
      label = current_statistic$raw_p_label[1],
      fontface = "bold",
      size = 3.6
    ) +
    scale_fill_manual(values = region_colors, drop = FALSE) +
    scale_x_discrete(
      labels = c(
        "Non-response region" = "Non-response",
        "Response region" = "Response"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.18))
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Standardized module score"
    ) +
    theme_manuscript +
    theme(
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 4),
      axis.text.x = element_text(
        angle = 12,
        hjust = 0.75,
        face = "bold",
        size = 8.2
      )
    )
}

make_figure3_gene_effect_plot <- function(
  gene_group_name,
  plot_title,
  show_legend = TRUE
) {
  current_data <- figure3_selected_gene_effects[
    as.character(gene_group) == gene_group_name
  ]

  if (nrow(current_data) == 0L) {
    stop("No gene-effect data are available for ", gene_group_name, ".")
  }

  current_plot <- ggplot(
    current_data,
    aes(
      x = log2FC_response_vs_nonresponse,
      y = gene_plot,
      fill = enriched_region
    )
  ) +
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linewidth = 0.5,
      linetype = 2
    ) +
    geom_col(width = 0.72, alpha = 0.92) +
    geom_text(
      aes(
        label = significance,
        hjust = ifelse(
          log2FC_response_vs_nonresponse >= 0,
          -0.18,
          1.18
        )
      ),
      fontface = "bold",
      size = 3.0
    ) +
    scale_fill_manual(values = region_colors, drop = FALSE) +
    scale_x_continuous(
      expand = expansion(mult = c(0.15, 0.15))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = plot_title,
      x = "log2 fold change (response / non-response)",
      y = NULL,
      fill = NULL
    ) +
    theme_manuscript +
    theme(
      legend.position = "top",
      axis.text.y = element_text(face = "bold", size = 7.8)
    )

  if (!isTRUE(show_legend)) {
    current_plot <- current_plot + theme(legend.position = "none")
  }

  current_plot
}

make_cancer_spatial_plot <- function(
  program_name,
  region_name,
  plot_title,
  show_legend = TRUE
) {
  current_background <- spatial_data[
    as.character(Region) == region_name
  ]

  current_cancer <- copy(
    cancer_score_data[
      program == program_name &
        as.character(Region) == region_name,
      .(
        cell,
        x_um,
        y_um,
        score_z,
        functional_population = "Cancer"
      )
    ]
  )

  current_niche <- copy(
    functional_data[
      program %in% c(
        "CAF_ECM_remodeling",
        "Myeloid_inflammatory"
      ) &
        as.character(Region) == region_name &
        distance_to_nearest_cancer_um <= niche_radius_um,
      .(
        cell,
        x_um,
        y_um,
        score_z,
        functional_population = fifelse(
          program == "CAF_ECM_remodeling",
          "ECM-CAF",
          "Inflammatory myeloid"
        )
      )
    ]
  )

  current_target <- rbindlist(
    list(
      current_cancer,
      current_niche
    ),
    use.names = TRUE,
    fill = TRUE
  )

  current_target[
    ,
    functional_population := factor(
      functional_population,
      levels = c(
        "Cancer",
        "ECM-CAF",
        "Inflammatory myeloid"
      )
    )
  ]

  current_background_grey <- current_background[
    !cell %in% current_target$cell
  ]

  common_limit <- score_limit(
    c(
      cancer_score_data[
        program == program_name,
        score_z
      ],
      functional_data[
        program %in% c(
          "CAF_ECM_remodeling",
          "Myeloid_inflammatory"
        ) &
          distance_to_nearest_cancer_um <= niche_radius_um,
        score_z
      ]
    )
  )

  current_plot <- ggplot() +
    geom_point(
      data = current_background_grey,
      aes(x = x_um, y = y_um),
      colour = "grey88",
      size = 0.30,
      alpha = 0.38
    ) +
    geom_point(
      data = current_niche,
      aes(
        x = x_um,
        y = y_um,
        fill = score_z,
        shape = functional_population,
        colour = functional_population
      ),
      stroke = 0.20,
      size = 1.52,
      alpha = 0.98
    ) +
    ## Draw cancer cells last so they remain visible beneath dense CAF/myeloid
    ## points. A white clearance halo and a black circular outline distinguish
    ## cancer cells without converting their score colour into a fixed colour.
    geom_point(
      data = current_cancer,
      aes(x = x_um, y = y_um),
      shape = 21,
      fill = "white",
      colour = "white",
      stroke = 0,
      size = 2.10,
      alpha = 0.98,
      show.legend = FALSE
    ) +
    geom_point(
      data = current_cancer,
      aes(x = x_um, y = y_um),
      shape = 21,
      fill = NA,
      colour = "#111111",
      stroke = 0.58,
      size = 1.78,
      alpha = 1,
      show.legend = FALSE
    ) +
    geom_point(
      data = current_cancer,
      aes(
        x = x_um,
        y = y_um,
        fill = score_z,
        shape = functional_population,
        colour = functional_population
      ),
      stroke = 0.12,
      size = 1.12,
      alpha = 0.98
    ) +
    scale_fill_gradient2(
      low = score_low_color,
      mid = "white",
      high = score_high_color,
      midpoint = 0,
      limits = c(-common_limit, common_limit),
      oob = squish,
      name = "Cell-type-specific\nstandardized score"
    ) +
    scale_colour_manual(
      values = c(
        "Cancer" = unname(celltype_colors[["Cancer"]]),
        "ECM-CAF" = unname(celltype_colors[["CAF_ECM"]]),
        "Inflammatory myeloid" =
          unname(celltype_colors[["Imm_Myeloid"]])
      ),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        "Cancer" = 21,
        "ECM-CAF" = 24,
        "Inflammatory myeloid" = 22
      ),
      breaks = c(
        "Cancer",
        "ECM-CAF",
        "Inflammatory myeloid"
      ),
      name = NULL
    ) +
    coord_equal(clip = "off") +
    labs(
      title = plot_title,
      x = "X coordinate (µm)",
      y = "Y coordinate (µm)",
      caption = paste0(
        "ECM-CAF and myeloid cells within ",
        niche_radius_um,
        " µm of cancer cells"
      )
    ) +
    theme_manuscript +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4),
      plot.caption = element_text(
        hjust = 1,
        size = 6.8,
        colour = "grey30"
      )
    ) +
    guides(
      fill = guide_colourbar(
        direction = "vertical",
        title.position = "top",
        barwidth = grid::unit(0.34, "cm"),
        barheight = grid::unit(3.0, "cm"),
        order = 1
      ),
      shape = guide_legend(
        direction = "vertical",
        ncol = 1,
        order = 2,
        override.aes = list(
          size = 3.2,
          colour = c(
            unname(celltype_colors[["Cancer"]]),
            unname(celltype_colors[["CAF_ECM"]]),
            unname(celltype_colors[["Imm_Myeloid"]])
          ),
          stroke = 0.45,
          alpha = 1,
          fill = "white"
        )
      )
    )

  if (!isTRUE(show_legend)) {
    current_plot <- current_plot + theme(legend.position = "none")
  } else {
    current_plot <- current_plot +
      theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.text = element_text(size = 8.8),
        legend.title = element_text(
          size = 8.8,
          face = "bold"
        ),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 3)
      )
  }

  current_plot
}

make_cancer_violin <- function(
  program_name,
  plot_title
) {
  current_data <- cancer_score_data[
    program == program_name
  ]

  current_statistic <- cancer_score_statistics[
    program == program_name
  ]

  if (nrow(current_data) == 0L || nrow(current_statistic) == 0L) {
    stop("No cancer-program data are available for ", program_name, ".")
  }

  y_maximum <- max(current_data$score_z, na.rm = TRUE)
  y_range <- diff(range(current_data$score_z, na.rm = TRUE))

  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 1
  }

  ggplot(
    current_data,
    aes(x = Region, y = score_z, fill = Region)
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      alpha = 0.55,
      linewidth = 0.4
    ) +
    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.45
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,
      size = 2.8,
      fill = "white",
      colour = "black"
    ) +
    annotate(
      geom = "text",
      x = 1.5,
      y = y_maximum + 0.12 * y_range,
      label = current_statistic$raw_p_label[1],
      fontface = "bold",
      size = 3.6
    ) +
    scale_fill_manual(values = region_colors, drop = FALSE) +
    scale_x_discrete(
      labels = c(
        "Non-response region" = "Non-response",
        "Response region" = "Response"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.18))
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Standardized module score"
    ) +
    theme_manuscript +
    theme(
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 2),
      axis.title.y = element_text(
        margin = margin(r = 3)
      ),
      axis.text.x = element_text(
        angle = 12,
        hjust = 0.75,
        face = "bold",
        size = 8.2
      )
    )
}

## ------------------------------------------------------------
## 14. Figure 3: ECM-CAF-myeloid spatial and functional states
## ------------------------------------------------------------

run_figure3_outputs <- !isTRUE(
  getOption(
    "crc.cosmx.figure4.only",
    FALSE
  )
)

if (run_figure3_outputs) {

figure3_required_programs <- c(
  "CAF_ECM_remodeling",
  "Myeloid_inflammatory",
  "Tcell_effector_activation",
  "Cytotoxic_NKT_cytotoxicity"
)

missing_figure3_programs <- setdiff(
  figure3_required_programs,
  unique(adjacent_functional_data$program)
)

if (length(missing_figure3_programs) > 0L) {
  stop(
    "Missing Figure 3 programs: ",
    paste(missing_figure3_programs, collapse = ", ")
  )
}

p3A <- make_combined_tme_spatial_plot(
  region_name = "Non-response region",
  plot_title = "Cancer–CAF–myeloid niche: Non-response",
  show_legend = FALSE,
  maximum_distance_um = niche_radius_um
)

p3B <- make_combined_tme_spatial_plot(
  region_name = "Response region",
  plot_title = "Cancer–CAF–myeloid niche: Response",
  show_legend = TRUE,
  maximum_distance_um = niche_radius_um
)

p3C <- make_function_violin(
  program_name = "CAF_ECM_remodeling",
  plot_title = paste0(
    "ECM-CAF remodeling (≤",
    primary_radius_um,
    " µm)"
  )
)

p3D <- make_function_violin(
  program_name = "Myeloid_inflammatory",
  plot_title = paste0(
    "Myeloid inflammation (≤",
    primary_radius_um,
    " µm)"
  )
)

p3E <- make_combined_immune_spatial_plot(
  region_name = "Non-response region",
  plot_title = "Tumor-adjacent immune states: Non-response",
  show_legend = FALSE
)

p3F <- make_combined_immune_spatial_plot(
  region_name = "Response region",
  plot_title = "Tumor-adjacent immune states: Response",
  show_legend = TRUE
)

p3G <- make_function_violin(
  program_name = "Tcell_effector_activation",
  plot_title = "T-cell effector/activation"
)

p3H <- make_function_violin(
  program_name = "Cytotoxic_NKT_cytotoxicity",
  plot_title = "Cytotoxic NK/T activity"
)

p3I <- make_figure3_gene_effect_plot(
  gene_group_name = "ECM-CAF remodeling",
  plot_title = "ECM-CAF remodeling genes",
  show_legend = TRUE
)

p3J <- make_figure3_gene_effect_plot(
  gene_group_name = "Myeloid inflammation",
  plot_title = "Myeloid inflammatory genes",
  show_legend = FALSE
)

p3K <- make_figure3_gene_effect_plot(
  gene_group_name = "T-cell effector/activation",
  plot_title = "T-cell effector/activation genes",
  show_legend = FALSE
)

p3L <- make_figure3_gene_effect_plot(
  gene_group_name = "Cytotoxic NK/T activity",
  plot_title = "Cytotoxic NK/T genes",
  show_legend = FALSE
)

## Final Figure 3 assembly is performed after the distance analysis below.

## ------------------------------------------------------------
## 15. Distance-dependent CAF-myeloid niche and final Figure 3
## ------------------------------------------------------------

distance_programs <- c(
  "CAF_ECM_remodeling",
  "Myeloid_inflammatory"
)

missing_distance_programs <- setdiff(
  distance_programs,
  unique(functional_data$program)
)

if (length(missing_distance_programs) > 0L) {
  stop(
    "Missing distance-analysis programs: ",
    paste(missing_distance_programs, collapse = ", ")
  )
}

distance_bin_upper_values <- primary_radius_um * seq_len(4L)

distance_bin_breaks <- c(
  0,
  distance_bin_upper_values,
  Inf
)

distance_bin_labels <- c(
  paste0(
    "0–",
    formatC(
      distance_bin_upper_values[1],
      format = "f",
      digits = 2
    )
  ),
  paste0(
    formatC(
      distance_bin_upper_values[1:3],
      format = "f",
      digits = 2
    ),
    "–",
    formatC(
      distance_bin_upper_values[2:4],
      format = "f",
      digits = 2
    )
  ),
  paste0(
    ">",
    formatC(
      distance_bin_upper_values[4],
      format = "f",
      digits = 2
    )
  )
)

distance_abundance_base <- copy(noncancer_distance_data)

distance_abundance_base[
  ,
  distance_bin := cut(
    distance_to_nearest_cancer_um,
    breaks = distance_bin_breaks,
    labels = distance_bin_labels,
    include.lowest = TRUE,
    right = FALSE
  )
]

distance_abundance_base <- distance_abundance_base[
  !is.na(distance_bin)
]

distance_denominator <- distance_abundance_base[
  ,
  .(
    all_noncancer_cells = as.integer(.N)
  ),
  by = .(
    Region,
    distance_bin
  )
]

distance_population_definitions <- data.table(
  population = c(
    "ECM-CAF",
    "Inflammatory myeloid"
  ),
  target_celltype = c(
    "CAF_ECM",
    "Imm_Myeloid"
  )
)

distance_population_count_list <- vector(
  mode = "list",
  length = nrow(distance_population_definitions)
)

for (
  population_index in seq_len(
    nrow(distance_population_definitions)
  )
) {
  current_definition <- distance_population_definitions[
    population_index
  ]

  current_counts <- distance_abundance_base[
    Final_celltype == current_definition$target_celltype,
    .(
      target_cells = as.integer(.N)
    ),
    by = .(
      Region,
      distance_bin
    )
  ]

  current_counts[
    ,
    population := current_definition$population
  ]

  distance_population_count_list[[population_index]] <-
    current_counts
}

distance_population_counts <- rbindlist(
  distance_population_count_list,
  use.names = TRUE,
  fill = TRUE
)

distance_abundance_grid <- CJ(
  Region = region_levels,
  distance_bin = distance_bin_labels,
  population = distance_population_definitions$population,
  unique = TRUE
)

distance_abundance_summary <- merge(
  distance_abundance_grid,
  distance_denominator[
    ,
    .(
      Region = as.character(Region),
      distance_bin = as.character(distance_bin),
      all_noncancer_cells
    )
  ],
  by = c(
    "Region",
    "distance_bin"
  ),
  all.x = TRUE,
  sort = FALSE
)

distance_abundance_summary <- merge(
  distance_abundance_summary,
  distance_population_counts[
    ,
    .(
      Region = as.character(Region),
      distance_bin = as.character(distance_bin),
      population,
      target_cells
    )
  ],
  by = c(
    "Region",
    "distance_bin",
    "population"
  ),
  all.x = TRUE,
  sort = FALSE
)

distance_abundance_summary[
  is.na(all_noncancer_cells),
  all_noncancer_cells := 0L
]

distance_abundance_summary[
  is.na(target_cells),
  target_cells := 0L
]

distance_abundance_summary[
  ,
  proportion_percent := fifelse(
    all_noncancer_cells > 0L,
    100 * as.numeric(target_cells) /
      as.numeric(all_noncancer_cells),
    NA_real_
  )
]

distance_abundance_summary[
  ,
  `:=`(
    Region = factor(
      Region,
      levels = region_levels
    ),
    distance_bin = factor(
      distance_bin,
      levels = distance_bin_labels
    )
  )
]

functional_distance_data <- copy(
  functional_data[
    program %in% distance_programs
  ]
)

functional_distance_data[
  ,
  distance_bin := cut(
    distance_to_nearest_cancer_um,
    breaks = distance_bin_breaks,
    labels = distance_bin_labels,
    include.lowest = TRUE,
    right = FALSE
  )
]

functional_distance_bin_summary <- functional_distance_data[
  !is.na(distance_bin),
  .(
    cells = as.integer(.N),
    mean_score_z = as.numeric(
      mean(score_z, na.rm = TRUE)
    ),
    median_score_z = as.numeric(
      median(score_z, na.rm = TRUE)
    ),
    sd_score_z = as.numeric(
      sd(score_z, na.rm = TRUE)
    ),
    se_score_z = as.numeric(
      sd(score_z, na.rm = TRUE) / sqrt(.N)
    )
  ),
  by = .(
    program,
    program_display,
    Region,
    distance_bin
  )
]

distance_score_statistics <- functional_distance_data[
  ,
  {
    current_test <- safe_spearman(
      distance_to_nearest_cancer_um,
      score_z
    )

    .(
      cells = as.integer(current_test$n),
      spearman_rho = as.numeric(current_test$rho),
      raw_p = as.numeric(current_test$raw_p)
    )
  },
  by = .(
    program,
    program_display,
    Region
  )
]

distance_score_statistics[
  ,
  FDR := p.adjust(
    raw_p,
    method = "BH"
  )
]

distance_score_statistics[
  ,
  correlation_label := paste0(
    "rho = ",
    fifelse(
      is.finite(spearman_rho),
      formatC(
        spearman_rho,
        format = "f",
        digits = 2
      ),
      "NA"
    ),
    "\n",
    vapply(
      raw_p,
      format_raw_p,
      character(1)
    )
  )
]

distance_heatmap_data <- rbindlist(
  list(
    distance_abundance_summary[
      population == "ECM-CAF",
      .(
        Region,
        distance_bin,
        feature = "ECM-CAF abundance",
        raw_value = proportion_percent
      )
    ],
    functional_distance_bin_summary[
      program == "CAF_ECM_remodeling",
      .(
        Region,
        distance_bin,
        feature = "ECM-CAF remodeling score",
        raw_value = mean_score_z
      )
    ],
    distance_abundance_summary[
      population == "Inflammatory myeloid",
      .(
        Region,
        distance_bin,
        feature = "Myeloid-cell abundance",
        raw_value = proportion_percent
      )
    ],
    functional_distance_bin_summary[
      program == "Myeloid_inflammatory",
      .(
        Region,
        distance_bin,
        feature = "Myeloid inflammatory score",
        raw_value = mean_score_z
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

distance_heatmap_data[
  ,
  scaled_value := standardize_vector(raw_value),
  by = feature
]

distance_heatmap_data[
  ,
  feature := factor(
    feature,
    levels = c(
      "ECM-CAF abundance",
      "ECM-CAF remodeling score",
      "Myeloid-cell abundance",
      "Myeloid inflammatory score"
    )
  )
]

distance_heatmap_data[
  ,
  `:=`(
    Region = factor(
      as.character(Region),
      levels = region_levels
    ),
    distance_bin = factor(
      as.character(distance_bin),
      levels = distance_bin_labels
    )
  )
]

fwrite(
  distance_abundance_summary,
  file.path(
    table3_dir,
    "Figure3_distance_bin_celltype_abundance.csv"
  )
)

fwrite(
  functional_distance_bin_summary,
  file.path(
    table3_dir,
    "Figure3_distance_bin_functional_scores.csv"
  )
)

fwrite(
  distance_score_statistics,
  file.path(
    table3_dir,
    "Figure3_distance_score_correlations.csv"
  )
)

fwrite(
  distance_heatmap_data,
  file.path(
    table3_dir,
    "Figure3_distance_heatmap_data.csv"
  )
)

distance_plot_max_um <- max(
  primary_radius_um * 4,
  as.numeric(
    quantile(
      functional_distance_data$distance_to_nearest_cancer_um,
      probs = 0.99,
      na.rm = TRUE,
      names = FALSE
    )
  )
)

make_distance_abundance_plot <- function(
  population_name,
  plot_title,
  y_axis_label
) {
  current_data <- distance_abundance_summary[
    population == population_name
  ]

  ggplot(
    current_data,
    aes(
      x = distance_bin,
      y = proportion_percent,
      colour = Region,
      group = Region
    )
  ) +
    geom_line(
      linewidth = 1.0,
      na.rm = TRUE
    ) +
    geom_point(
      size = 2.8,
      stroke = 0.3,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = region_colors,
      drop = FALSE
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.12))
    ) +
    labs(
      title = plot_title,
      x = "Distance to nearest cancer cell (µm)",
      y = y_axis_label,
      colour = NULL
    ) +
    theme_manuscript +
    theme(
      legend.position = "top",
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.major.y = element_line(
        colour = "grey92",
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank()
    )
}

make_distance_score_plot <- function(
  program_name,
  plot_title,
  y_axis_label
) {
  current_data <- functional_distance_data[
    program == program_name &
      distance_to_nearest_cancer_um <= distance_plot_max_um
  ]

  current_statistics <- current_data[
    ,
    {
      current_test <- safe_spearman(
        distance_to_nearest_cancer_um,
        score_z
      )

      .(
        spearman_rho = as.numeric(current_test$rho),
        raw_p = as.numeric(current_test$raw_p)
      )
    },
    by = Region
  ]

  current_statistics[
    ,
    correlation_label := paste0(
      "rho = ",
      fifelse(
        is.finite(spearman_rho),
        formatC(
          spearman_rho,
          format = "f",
          digits = 2
        ),
        "NA"
      ),
      "\n",
      vapply(
        raw_p,
        format_raw_p,
        character(1)
      )
    )
  ]

  x_breaks <- sort(
    unique(
      c(
        0,
        primary_radius_um,
        pretty(
          c(0, distance_plot_max_um),
          n = 5
        )
      )
    )
  )

  x_breaks <- x_breaks[
    x_breaks >= 0 &
      x_breaks <= distance_plot_max_um
  ]

  ggplot(
    current_data,
    aes(
      x = distance_to_nearest_cancer_um,
      y = score_z,
      colour = Region
    )
  ) +
    geom_vline(
      xintercept = primary_radius_um,
      colour = "black",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_point(
      size = 0.48,
      alpha = 0.16,
      stroke = 0
    ) +
    geom_smooth(
      aes(fill = Region),
      method = "loess",
      formula = y ~ x,
      span = 0.80,
      se = TRUE,
      linewidth = 1.0,
      alpha = 0.15
    ) +
    geom_text(
      data = current_statistics,
      aes(
        x = Inf,
        y = Inf,
        label = correlation_label
      ),
      inherit.aes = FALSE,
      hjust = 1.08,
      vjust = 1.15,
      fontface = "bold",
      size = 3.1
    ) +
    scale_colour_manual(
      values = region_colors,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = region_colors,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      limits = c(0, distance_plot_max_um),
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    facet_wrap(
      facets = vars(Region),
      nrow = 1,
      drop = FALSE
    ) +
    labs(
      title = plot_title,
      x = "Distance to nearest cancer cell (µm)",
      y = y_axis_label,
      caption = paste0(
        "Dashed line: ",
        primary_radius_um,
        " µm tumor-adjacent boundary"
      )
    ) +
    guides(
      colour = "none",
      fill = "none"
    ) +
    theme_manuscript +
    theme(
      panel.grid.major = element_line(
        colour = "grey92",
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(
        hjust = 1,
        size = 7.5,
        colour = "grey30"
      )
    )
}

pS3A <- make_distance_abundance_plot(
  population_name = "ECM-CAF",
  plot_title = "ECM-CAF abundance by tumor distance",
  y_axis_label = "ECM-CAF among non-cancer cells (%)"
)

pS3B <- make_distance_score_plot(
  program_name = "CAF_ECM_remodeling",
  plot_title = "ECM-CAF remodeling across tumor distance",
  y_axis_label = "ECM-CAF remodeling score (z)"
)

pS3C <- make_distance_abundance_plot(
  population_name = "Inflammatory myeloid",
  plot_title = "Myeloid-cell abundance by tumor distance",
  y_axis_label = "Myeloid cells among non-cancer cells (%)"
)

pS3D <- make_distance_score_plot(
  program_name = "Myeloid_inflammatory",
  plot_title = "Myeloid inflammation across tumor distance",
  y_axis_label = "Myeloid inflammatory score (z)"
)

heatmap_limit <- score_limit(
  distance_heatmap_data$scaled_value
)

p3K_distance_heatmap <- ggplot(
  distance_heatmap_data,
  aes(
    x = distance_bin,
    y = feature,
    fill = scaled_value
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = formatC(
        raw_value,
        format = "f",
        digits = 2
      )
    ),
    size = 3.0,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = score_low_color,
    mid = "white",
    high = score_high_color,
    midpoint = 0,
    limits = c(-heatmap_limit, heatmap_limit),
    oob = squish,
    name = "Within-feature\nz-score"
  ) +
  facet_wrap(
    facets = vars(Region),
    nrow = 1,
    drop = FALSE
  ) +
  labs(
    title = "Coordinated CAF-myeloid distance profile",
    x = "Distance to nearest cancer cell (µm)",
    y = NULL
  ) +
  theme_manuscript +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    panel.grid = element_blank()
  )

## Distance plots used in the final main figure.
p3E_distance <- pS3A
p3F_distance <- pS3B
p3G_distance <- pS3C
p3H_distance <- pS3D

figure3_design <- paste(
  "AABBCCDD",
  "EEFFGGHH",
  "IIJJKKKK",
  sep = "\n"
)

figure3 <- wrap_plots(
  p3A,
  p3B,
  p3C,
  p3D,
  p3E_distance,
  p3F_distance,
  p3G_distance,
  p3H_distance,
  p3I,
  p3J,
  p3K_distance_heatmap,
  design = figure3_design,
  widths = rep(1, 8),
  heights = c(1, 1, 1.08),
  guides = "keep"
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.005, 0.995)
  )

print(figure3)

save_plot_formats(
  plot_object = figure3,
  output_directory = figure3_dir,
  base_filename = "Figure3",
  width = 22.5,
  height = 18.0,
  save_tiff = TRUE
)

figure3_panels <- list(
  p3A,
  p3B,
  p3C,
  p3D,
  p3E_distance,
  p3F_distance,
  p3G_distance,
  p3H_distance,
  p3I,
  p3J,
  p3K_distance_heatmap
)

save_panel_collection(
  panel_list = figure3_panels,
  output_directory = figure3_dir,
  prefix = "Figure3",
  width = c(
    8.2, 8.2, 7.0, 7.0,
    8.0, 12.0, 8.0, 12.0,
    9.0, 9.0, 14.0
  ),
  height = c(
    7.2, 7.2, 6.2, 6.2,
    6.4, 6.4, 6.4, 6.4,
    7.2, 7.2, 7.0
  )
)

## ------------------------------------------------------------
## 15B. Supplementary Figure 3: T-cell and cytotoxic NK/T states
## ------------------------------------------------------------

supplementary_figure3_design <- paste(
  "ABCD",
  "EEFF",
  sep = "\n"
)

supplementary_figure3 <- wrap_plots(
  p3E,
  p3F,
  p3G,
  p3H,
  p3K,
  p3L,
  design = supplementary_figure3_design,
  heights = c(0.92, 0.98),
  guides = "keep"
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.005, 0.995)
  )

print(supplementary_figure3)

save_plot_formats(
  plot_object = supplementary_figure3,
  output_directory = supplementary3_dir,
  base_filename = "Supplementary Figure3",
  width = 22.5,
  height = 11.2,
  save_tiff = TRUE
)

supplementary3_panels <- list(
  p3E,
  p3F,
  p3G,
  p3H,
  p3K,
  p3L
)

save_panel_collection(
  panel_list = supplementary3_panels,
  output_directory = supplementary3_dir,
  prefix = "Supplementary Figure3",
  width = c(8.0, 8.0, 7.0, 7.0, 9.0, 9.0),
  height = c(7.0, 7.0, 6.2, 6.2, 7.2, 7.2)
)

} else {
  cat(
    "\nFigure 4-only mode: skipping Figure 3 and ",
    "Supplementary Figure 3 outputs.\n"
  )
}

## ------------------------------------------------------------
## 16. Figure 4: residual cancer-cell programs
## ------------------------------------------------------------

figure4_required_programs <- c(
  "Cancer_IFN_adaptive_resistance",
  "Cancer_antigen_presentation"
)

missing_figure4_programs <- setdiff(
  figure4_required_programs,
  unique(cancer_score_data$program)
)

if (length(missing_figure4_programs) > 0L) {
  stop(
    "Missing Figure 4 programs: ",
    paste(missing_figure4_programs, collapse = ", ")
  )
}

p4A <- make_cancer_spatial_plot(
  program_name = "Cancer_IFN_adaptive_resistance",
  region_name = "Non-response region",
  plot_title = "Cancer–CAF–myeloid states: Non-response",
  show_legend = FALSE
)

p4B <- make_cancer_spatial_plot(
  program_name = "Cancer_IFN_adaptive_resistance",
  region_name = "Response region",
  plot_title = "Cancer–CAF–myeloid states: Response",
  show_legend = TRUE
)

p4C <- make_cancer_violin(
  program_name = "Cancer_IFN_adaptive_resistance",
  plot_title = "Cancer IFN/adaptive resistance"
)

p4D <- make_cancer_violin(
  program_name = "Cancer_antigen_presentation",
  plot_title = "Cancer antigen presentation"
)

p4E <- ggplot(
  selected_gene_effects,
  aes(
    x = log2FC_response_vs_nonresponse,
    y = gene_plot,
    fill = enriched_region
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "grey45",
    linewidth = 0.5,
    linetype = 2
  ) +
  geom_col(width = 0.72, alpha = 0.92) +
  geom_text(
    aes(
      label = significance,
      hjust = ifelse(
        log2FC_response_vs_nonresponse >= 0,
        -0.18,
        1.18
      )
    ),
    fontface = "bold",
    size = 3.1
  ) +
  scale_fill_manual(values = region_colors, drop = FALSE) +
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  facet_grid(
    rows = vars(gene_group),
    scales = "free_y",
    space = "free_y"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Regional cancer-cell gene programs",
    x = "log2 fold change (response / non-response)",
    y = NULL,
    fill = NULL
  ) +
  theme_manuscript +
  theme(
    legend.position = "top",
    axis.text.y = element_text(face = "bold", size = 8.2),
    strip.text.y = element_text(angle = 0, size = 8.2)
  )

maximum_local_ecm_caf_count <- max(
  ecm_ifn_coupling_data$local_ecm_caf_count,
  na.rm = TRUE
)

candidate_count_breaks <- c(
  0, 1, 2, 5, 10, 20, 50, 100, 200
)

local_count_breaks <- candidate_count_breaks[
  candidate_count_breaks <= maximum_local_ecm_caf_count
]

if (
  maximum_local_ecm_caf_count > 0L &&
    length(local_count_breaks) < 2L
) {
  local_count_breaks <- c(
    0,
    maximum_local_ecm_caf_count
  )
}

ecm_activity_limit <- score_limit(
  ecm_ifn_coupling_data[
    local_ecm_caf_count > 0L,
    local_mean_ecm_score_z
  ]
)

pS4B_ecm_ifn <- ggplot(
  ecm_ifn_coupling_data,
  aes(
    x = log1p_local_ecm_caf_count,
    y = cancer_ifn_score_z
  )
) +
  geom_point(
    aes(
      fill = local_mean_ecm_score_z,
      colour = Region
    ),
    shape = 21,
    size = 1.25,
    stroke = 0.22,
    alpha = 0.48,
    position = position_jitter(
      width = 0.035,
      height = 0,
      seed = 20261008
    )
  ) +
  geom_smooth(
    aes(colour = Region),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.85,
    alpha = 0.13
  ) +
  geom_text(
    data = ecm_ifn_correlation_statistics,
    aes(
      x = Inf,
      y = Inf,
      label = count_correlation_label
    ),
    inherit.aes = FALSE,
    hjust = 1.08,
    vjust = 1.15,
    fontface = "bold",
    size = 3.0
  ) +
  scale_x_continuous(
    breaks = log1p(local_count_breaks),
    labels = local_count_breaks,
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  scale_fill_gradient2(
    low = score_low_color,
    mid = "white",
    high = score_high_color,
    midpoint = 0,
    limits = c(-ecm_activity_limit, ecm_activity_limit),
    oob = squish,
    na.value = "grey82",
    name = paste0(
      "Mean local ECM-CAF\n",
      "remodeling score (z)"
    )
  ) +
  scale_colour_manual(values = region_colors, drop = FALSE) +
  facet_wrap(
    facets = vars(Region),
    nrow = 1,
    drop = FALSE
  ) +
  labs(
    title = "Cancer IFN signaling versus local ECM-CAF burden",
    x = paste0(
      "ECM-CAF cells within ",
      primary_radius_um,
      " µm"
    ),
    y = "Cancer IFN/adaptive-resistance score (z)"
  ) +
  guides(colour = "none") +
  theme_manuscript +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(
      colour = "grey92",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8.2)
  )

## Selected interactions linking the Figure 3 CAF-myeloid niche to
## the regional cancer-cell programs in Figure 4.
main_interaction_priority <- c(
  "Cancer__CSF1__CSF1R__Imm_Myeloid",
  "Cancer__CD47__SIRPA__Imm_Myeloid",
  "CAF_ECM__FN1__CD44__Cancer",
  "CAF_ECM__TGFB1__TGFBR2__Cancer",
  "Imm_Myeloid__IL1B__IL1R1__Cancer",
  "Imm_Myeloid__SPP1__CD44__Cancer",
  "Cancer__CXCL9__CXCR3__T/NK",
  "Cancer__CXCL10__CXCR3__T/NK"
)

main_interaction_plot_data <- copy(
  interaction_statistics[
    interaction_id %in% main_interaction_priority &
      eligible_for_main_figure == TRUE
  ]
)

main_interaction_plot_data[
  ,
  main_priority := match(
    interaction_id,
    main_interaction_priority
  )
]

minimum_main_interactions <- min(
  6L,
  nrow(interaction_plot_data)
)

if (nrow(main_interaction_plot_data) < minimum_main_interactions) {
  fallback_interactions <- copy(
    interaction_plot_data[
      !interaction_id %in%
        main_interaction_plot_data$interaction_id
    ]
  )

  fallback_interactions[
    ,
    main_priority :=
      length(main_interaction_priority) + seq_len(.N)
  ]

  main_interaction_plot_data <- rbindlist(
    list(
      main_interaction_plot_data,
      fallback_interactions
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

setorder(
  main_interaction_plot_data,
  main_priority,
  priority
)

main_interaction_plot_data <- unique(
  main_interaction_plot_data,
  by = "interaction_id"
)

main_interaction_plot_data <- head(
  main_interaction_plot_data,
  8L
)

if (nrow(main_interaction_plot_data) == 0L) {
  stop(
    "No cancer-centered spatial interaction was available for Figure 4F."
  )
}

main_interaction_plot_data[
  ,
  interaction_display := paste0(
    pair_label,
    "  (",
    sender_group,
    "→",
    receiver_group,
    ")"
  )
]

main_interaction_plot_data[
  ,
  interaction_plot := factor(
    interaction_display,
    levels = rev(interaction_display)
  )
]

main_interaction_codetection_long <- rbindlist(
  list(
    main_interaction_plot_data[
      ,
      .(
        interaction_id,
        interaction_plot,
        Region = "Non-response region",
        codetection_percent =
          100 * codetection_nonresponse
      )
    ],
    main_interaction_plot_data[
      ,
      .(
        interaction_id,
        interaction_plot,
        Region = "Response region",
        codetection_percent =
          100 * codetection_response
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

main_interaction_codetection_long[
  ,
  Region := factor(Region, levels = region_levels)
]

fwrite(
  main_interaction_plot_data,
  file.path(
    table4_dir,
    "Figure4_selected_spatial_interactions.csv"
  )
)

p4F <- ggplot(
  main_interaction_codetection_long,
  aes(
    x = codetection_percent,
    y = interaction_plot,
    colour = Region
  )
) +
  geom_segment(
    data = main_interaction_plot_data,
    aes(
      x = 100 * codetection_nonresponse,
      xend = 100 * codetection_response,
      y = interaction_plot,
      yend = interaction_plot
    ),
    inherit.aes = FALSE,
    colour = "grey68",
    linewidth = 0.85
  ) +
  geom_point(
    size = 4.0,
    alpha = 0.96
  ) +
  scale_colour_manual(
    values = region_colors,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.12))
  ) +
  labs(
    title = "Cancer-centered spatial signaling niches",
    x = paste0(
      "Ligand-positive sender cells with a receptor-positive neighbor within ",
      primary_radius_um,
      " µm (%)"
    ),
    y = NULL,
    colour = NULL
  ) +
  theme_manuscript +
  theme(
    legend.position = "top",
    axis.text.y = element_text(
      face = "bold",
      size = 8.0
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## Build the two rows independently. This prevents the long y-axis labels in
## panel F from forcing an oversized left alignment space into panel C.
figure4_tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 18),
  plot.tag.position = c(0.005, 0.995)
)

figure4_top_row <- wrap_plots(
  p4A + labs(tag = "A"),
  p4B + labs(tag = "B"),
  p4C + labs(tag = "C"),
  p4D + labs(tag = "D"),
  nrow = 1,
  widths = c(1.00, 1.27, 0.88, 0.88),
  guides = "keep"
) &
  figure4_tag_theme

figure4_bottom_row <- wrap_plots(
  p4E + labs(tag = "E"),
  p4F + labs(tag = "F"),
  nrow = 1,
  widths = c(1, 1),
  guides = "keep"
) &
  figure4_tag_theme

figure4 <- wrap_plots(
  figure4_top_row,
  figure4_bottom_row,
  ncol = 1,
  heights = c(1, 1.18)
)

print(figure4)

save_plot_formats(
  plot_object = figure4,
  output_directory = figure4_dir,
  base_filename = "Figure4",
  width = 22.0,
  height = 13.6,
  save_tiff = TRUE
)

figure4_panels <- list(p4A, p4B, p4C, p4D, p4E, p4F)

save_panel_collection(
  panel_list = figure4_panels,
  output_directory = figure4_dir,
  prefix = "Figure4",
  width = c(8.0, 7.6, 7.6, 7.6, 8.2, 9.2),
  height = c(7.0, 7.0, 6.2, 6.2, 8.2, 7.2)
)

## ------------------------------------------------------------
## 17. Supplementary Figure 4: extended niche-coupling analyses
## ------------------------------------------------------------

interaction_codetection_long <- rbindlist(
  list(
    interaction_plot_data[
      ,
      .(
        interaction_id,
        interaction_class,
        pair_label,
        priority,
        Region = "Non-response region",
        codetection_percent = 100 * codetection_nonresponse
      )
    ],
    interaction_plot_data[
      ,
      .(
        interaction_id,
        interaction_class,
        pair_label,
        priority,
        Region = "Response region",
        codetection_percent = 100 * codetection_response
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

interaction_codetection_long[
  ,
  Region := factor(Region, levels = region_levels)
]
interaction_codetection_long[
  ,
  interaction_class := factor(
    interaction_class,
    levels = unique(interaction_candidates$interaction_class)
  )
]
interaction_codetection_long[
  ,
  pair_plot := factor(
    pair_label,
    levels = levels(interaction_plot_data$pair_plot)
  )
]

pS4A <- ggplot(
  interaction_codetection_long,
  aes(
    x = codetection_percent,
    y = pair_plot,
    colour = Region
  )
) +
  geom_segment(
    data = interaction_plot_data,
    aes(
      x = 100 * codetection_nonresponse,
      xend = 100 * codetection_response,
      y = pair_plot,
      yend = pair_plot
    ),
    inherit.aes = FALSE,
    colour = "grey65",
    linewidth = 0.75
  ) +
  geom_point(
    size = 4.1,
    alpha = 0.95
  ) +
  facet_grid(
    rows = vars(interaction_class),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_colour_manual(values = region_colors, drop = FALSE) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  labs(
    title = "Spatially co-detected cancer-centered ligand-receptor pairs",
    x = paste0(
      "Sender cells with ligand and a receptor-positive neighbor within ",
      primary_radius_um,
      " µm (%)"
    ),
    y = NULL,
    colour = NULL
  ) +
  theme_manuscript +
  theme(
    legend.position = "top",
    axis.text.y = element_text(face = "bold", size = 8.5),
    strip.text.y = element_text(
      face = "bold",
      angle = 0,
      size = 8.2
    )
  )

supplementary_figure4 <- pS4A / pS4B_ecm_ifn +
  plot_layout(
    heights = c(1.08, 0.92),
    guides = "keep"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.005, 0.995)
  )

print(supplementary_figure4)

save_plot_formats(
  plot_object = supplementary_figure4,
  output_directory = supplementary4_dir,
  base_filename = "Supplementary Figure4",
  width = 17.0,
  height = 15.5,
  save_tiff = TRUE
)

save_panel_collection(
  panel_list = list(
    pS4A,
    pS4B_ecm_ifn
  ),
  output_directory = supplementary4_dir,
  prefix = "Supplementary Figure4",
  width = c(12.5, 12.5),
  height = c(10.0, 7.4)
)

## ------------------------------------------------------------
## 18. Save the scored object and reproducibility information
## ------------------------------------------------------------

if (isTRUE(save_scored_object)) {
  scored_object_file <- file.path(
    object_dir,
    "CRC_CosMx_postRT_Figure3_Figure4_scored_object.rds"
  )

  saveRDS(
    full_object,
    scored_object_file,
    compress = FALSE
  )
}

session_info_file <- file.path(
  table4_dir,
  "Figure3_Figure4_sessionInfo.txt"
)

writeLines(
  capture.output(sessionInfo()),
  con = session_info_file
)

cat("\n========================================\n")
cat("Figure 3 and Figure 4 analysis completed\n")
cat("========================================\n")
cat("Post-RT sample:", analysis_sample, "\n")
cat("Non-response region: FOV4\n")
cat("Response region: FOV6\n")
cat("Tumor-adjacent radius:", primary_radius_um, "µm\n")
cat("Extended CAF-myeloid niche radius:", niche_radius_um, "µm\n")
cat("\nFigure 3 primary tumor-adjacent CAF-myeloid statistics:\n")
print(
  functional_statistics[
    program %in% c(
      "CAF_ECM_remodeling",
      "Myeloid_inflammatory"
    ),
    .(
      program,
      Final_celltype,
      nonresponse_cells,
      response_cells,
      nonresponse_median_z,
      response_median_z,
      response_minus_nonresponse,
      raw_p,
      raw_p_label
    )
  ]
)
cat("\nFigure 3 distance-dependent CAF-myeloid correlations:\n")
print(
  distance_score_statistics[
    program %in% c(
      "CAF_ECM_remodeling",
      "Myeloid_inflammatory"
    ),
    .(
      program,
      Region,
      cells,
      spearman_rho,
      raw_p,
      FDR
    )
  ]
)
cat("\nSupplementary Figure 3 immune-state statistics:\n")
print(
  functional_statistics[
    program %in% c(
      "Tcell_effector_activation",
      "Cytotoxic_NKT_cytotoxicity"
    ),
    .(
      program,
      Final_celltype,
      nonresponse_cells,
      response_cells,
      nonresponse_median_z,
      response_median_z,
      response_minus_nonresponse,
      raw_p,
      raw_p_label
    )
  ]
)
cat("\nFigure 4F selected cancer-centered spatial interactions:\n")
print(
  main_interaction_plot_data[
    ,
    .(
      interaction_class,
      sender_group,
      pair_label,
      receiver_group,
      codetection_nonresponse,
      codetection_response,
      rank_biserial_response_vs_nonresponse,
      raw_p,
      FDR
    )
  ]
)
cat("\nSupplementary Figure 4B cancer IFN versus local ECM-CAF correlations:\n")
print(
  ecm_ifn_correlation_statistics[
    ,
    .(
      Region,
      cancer_cells,
      cancer_cells_with_local_ecm_caf,
      median_local_ecm_caf_count,
      rho_ifn_vs_ecm_caf_count,
      raw_p_ifn_vs_ecm_caf_count,
      rho_ifn_vs_mean_ecm_activity,
      raw_p_ifn_vs_mean_ecm_activity
    )
  ]
)
cat(
  "Ligand-receptor pairs available in panel:",
  nrow(available_interactions),
  "of",
  nrow(interaction_candidates),
  "\n"
)
cat(
  "Ligand-receptor pairs shown in Supplementary Figure 4:",
  nrow(interaction_plot_data),
  "\n"
)
cat("\nFigure 3:\n", file.path(figure3_dir, "Figure3.pdf"), "\n")
cat("\nSupplementary Figure 3:\n")
cat(file.path(supplementary3_dir, "Supplementary Figure3.pdf"), "\n")
cat("\nFigure 4:\n", file.path(figure4_dir, "Figure4.pdf"), "\n")
cat("\nSupplementary Figure 4:\n")
cat(file.path(supplementary4_dir, "Supplementary Figure4.pdf"), "\n")
cat("\nFigure 3 tables:\n", table3_dir, "\n")
cat("\nFigure 4 tables:\n", table4_dir, "\n")

if (isTRUE(save_scored_object)) {
  cat("\nScored object:\n", scored_object_file, "\n")
}

cat("========================================\n")

