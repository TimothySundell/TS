#' Find sample-level Cell Ranger output folders
#'
#' @description
#' Recursively searches a project folder for Cell Ranger matrix directories and
#' infers the corresponding sample names. The function is designed for common
#' `cellranger multi` layouts such as `outs/per_sample_outs/<sample>` and for
#' project folders organised as `project/<sample>/cellranger_outs`.
#'
#' The filtered gene expression matrix is detected at sample level. V(D)J
#' directories are detected at run level, for example `outs/vdj_b` and
#' `outs/vdj_t`, because `cellranger multi` stores the unfiltered
#' `all_contig*` files there rather than in `per_sample_outs/<sample>`.
#'
#' @param input_dir Character scalar. Project or Cell Ranger output directory to
#'   search.
#' @param sample_regex Optional regular expression used to extract sample names
#'   from matrix directory paths. If the regular expression contains a capture
#'   group, the first captured group is used; otherwise the full match is used.
#' @param matrix_dir_names Character vector of matrix directory names to search
#'   for. Defaults to Cell Ranger's `filtered_feature_bc_matrix` and
#'   `sample_filtered_feature_bc_matrix`.
#' @param vdj_dir_names Character vector of V(D)J directory names to detect.
#'
#' @return A tibble with one row per detected sample matrix and columns
#'   `sample_id`, `sample_root`, `outs_dir`, `matrix_dir`, `vdj_b_dir`,
#'   `vdj_t_dir`, and `detection_note`.
#'
#' @examples
#' project <- file.path(tempdir(), "cellranger_project")
#' matrix_dir <- file.path(project, "sampleA", "cellranger_outs", "filtered_feature_bc_matrix")
#' dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(matrix_dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")))
#'
#' TS_find_cellranger_samples(project)
#'
#' @export
TS_find_cellranger_samples <- function(
  input_dir,
  sample_regex = NULL,
  matrix_dir_names = c(
    "filtered_feature_bc_matrix",
    "sample_filtered_feature_bc_matrix"
  ),
  vdj_dir_names = c("vdj_b", "vdj_t")
) {
  .TS_check_character_scalar(input_dir, "input_dir")

  if (!dir.exists(input_dir)) {
    stop("input_dir does not exist: ", input_dir, call. = FALSE)
  }

  if (!is.null(sample_regex)) {
    .TS_check_character_scalar(sample_regex, "sample_regex")
  }

  if (!is.character(matrix_dir_names) || length(matrix_dir_names) == 0) {
    stop(
      "matrix_dir_names must be a non-empty character vector.",
      call. = FALSE
    )
  }

  if (!is.character(vdj_dir_names) || length(vdj_dir_names) == 0) {
    stop("vdj_dir_names must be a non-empty character vector.", call. = FALSE)
  }

  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  all_dirs <- .TS_list_dirs(input_dir)
  matrix_dirs <- all_dirs[basename(all_dirs) %in% matrix_dir_names]
  matrix_dirs <- matrix_dirs[!.TS_is_combined_multi_matrix_dir(matrix_dirs)]

  if (length(matrix_dirs) == 0) {
    return(tibble::tibble(
      sample_id = character(),
      sample_root = character(),
      outs_dir = character(),
      matrix_dir = character(),
      vdj_b_dir = character(),
      vdj_t_dir = character(),
      detection_note = character()
    ))
  }

  samples <- lapply(matrix_dirs, function(matrix_dir) {
    inferred <- .TS_infer_cellranger_sample(
      matrix_dir,
      sample_regex = sample_regex
    )
    vdj_dirs <- .TS_find_raw_vdj_dirs(
      inferred$outs_dir,
      vdj_dir_names = vdj_dir_names
    )

    tibble::tibble(
      sample_id = inferred$sample_id,
      sample_root = inferred$sample_root,
      outs_dir = inferred$outs_dir,
      matrix_dir = matrix_dir,
      vdj_b_dir = .TS_first_or_na(vdj_dirs[basename(vdj_dirs) == "vdj_b"]),
      vdj_t_dir = .TS_first_or_na(vdj_dirs[basename(vdj_dirs) == "vdj_t"]),
      detection_note = inferred$detection_note
    )
  })

  dplyr::bind_rows(samples)
}

#' Plan collection of selected Cell Ranger files
#'
#' @description
#' Builds a copy plan for selected Cell Ranger files without changing the file
#' system. Matrix files are planned from sample-level filtered matrix output
#' into per-sample `filtered_feature_bc_matrix` directories with their original
#' filenames, so they remain directly readable by Seurat. V(D)J files are
#' planned from run-level unfiltered `outs/vdj_b` and `outs/vdj_t` directories
#' into per-sample `vdj_b` and `vdj_t` directories with the sample name
#' prepended.
#'
#' @param input_dir Character scalar. Project or Cell Ranger output directory to
#'   search.
#' @param dest_dir Character scalar. Destination directory for the collected
#'   sample folders.
#' @param sample_regex Optional regular expression used to extract sample names.
#'   See [TS_find_cellranger_samples()].
#' @param overwrite Logical. If `FALSE`, the function errors when planned target
#'   files already exist.
#' @param strict Logical. If `TRUE`, existing raw V(D)J directories must contain
#'   all requested `all_contig*` files. If `FALSE`, missing V(D)J files are
#'   recorded in the plan and available files can still be copied.
#' @param ... Additional arguments passed to [TS_find_cellranger_samples()], such
#'   as `matrix_dir_names` or `vdj_dir_names`.
#'
#' @return A tibble with one row per planned or missing file. Important columns
#'   include `sample_id`, `file_group`, `source_path`, `dest_path`,
#'   `source_exists`, `copy`, and `status`.
#'
#' @examples
#' project <- file.path(tempdir(), "cellranger_project_plan")
#' matrix_dir <- file.path(project, "sampleA", "cellranger_outs", "filtered_feature_bc_matrix")
#' dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(matrix_dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")))
#'
#' plan <- TS_plan_cellranger_file_collection(
#'   input_dir = project,
#'   dest_dir = file.path(tempdir(), "cellranger_transfer_plan")
#' )
#' plan
#'
#' @export
TS_plan_cellranger_file_collection <- function(
  input_dir,
  dest_dir,
  sample_regex = NULL,
  overwrite = FALSE,
  strict = TRUE,
  ...
) {
  .TS_check_character_scalar(dest_dir, "dest_dir")
  .TS_check_logical_scalar(overwrite, "overwrite")
  .TS_check_logical_scalar(strict, "strict")

  samples <- TS_find_cellranger_samples(
    input_dir = input_dir,
    sample_regex = sample_regex,
    ...
  )

  if (nrow(samples) == 0) {
    stop(
      "No sample-level Cell Ranger matrix directories were found.",
      call. = FALSE
    )
  }

  .TS_validate_sample_ids(samples$sample_id)

  duplicated_samples <- unique(samples$sample_id[duplicated(samples$sample_id)])
  if (length(duplicated_samples) > 0) {
    stop(
      "Multiple matrix directories were detected for the same sample ID:\n",
      paste(duplicated_samples, collapse = "\n"),
      call. = FALSE
    )
  }

  .TS_validate_cellranger_vdj_sample_scope(samples)

  dest_dir <- normalizePath(dest_dir, winslash = "/", mustWork = FALSE)

  matrix_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  vdj_files <- c("all_contig_annotations.csv", "all_contig.fasta")

  plan <- dplyr::bind_rows(lapply(seq_len(nrow(samples)), function(i) {
    sample_row <- samples[i, , drop = FALSE]
    .TS_plan_one_cellranger_sample(
      sample_row = sample_row,
      dest_dir = dest_dir,
      matrix_files = matrix_files,
      vdj_files = vdj_files
    )
  }))

  missing_matrix <- dplyr::filter(
    plan,
    .data$required,
    !.data$source_exists
  )
  if (nrow(missing_matrix) > 0) {
    stop(
      "Some required filtered_feature_bc_matrix files are missing:\n",
      paste(
        paste0(missing_matrix$sample_id, ": ", missing_matrix$source_path),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  missing_vdj <- dplyr::filter(
    plan,
    !.data$required,
    !.data$source_exists
  )
  if (isTRUE(strict) && nrow(missing_vdj) > 0) {
    stop(
      "Some V(D)J directories are present but missing requested files:\n",
      paste(
        paste0(
          missing_vdj$sample_id,
          " ",
          missing_vdj$file_group,
          ": ",
          missing_vdj$source_path
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  dup_targets <- dplyr::filter(plan, .data$copy)
  dup_targets <- dplyr::count(dup_targets, .data$dest_path, name = "n")
  dup_targets <- dplyr::filter(dup_targets, .data$n > 1)

  if (nrow(dup_targets) > 0) {
    stop(
      "Multiple source files would be copied to the same target path:\n",
      paste(dup_targets$dest_path, collapse = "\n"),
      call. = FALSE
    )
  }

  same_paths <- dplyr::filter(
    plan,
    .data$copy,
    normalizePath(.data$source_path, winslash = "/", mustWork = FALSE) ==
      normalizePath(.data$dest_path, winslash = "/", mustWork = FALSE)
  )
  if (nrow(same_paths) > 0) {
    stop(
      "Some source and destination paths are identical:\n",
      paste(same_paths$source_path, collapse = "\n"),
      call. = FALSE
    )
  }

  existing_targets <- dplyr::filter(
    plan,
    .data$copy,
    file.exists(.data$dest_path)
  )
  if (nrow(existing_targets) > 0 && !isTRUE(overwrite)) {
    stop(
      "Target files already exist and overwrite = FALSE:\n",
      paste(existing_targets$dest_path, collapse = "\n"),
      call. = FALSE
    )
  }

  plan
}

#' Collect selected Cell Ranger files into per-sample folders
#'
#' @description
#' Copies selected Cell Ranger files into a transfer-friendly destination
#' directory. Each sample receives its own folder. Filtered matrix filenames are
#' kept unchanged for Seurat compatibility, and unfiltered V(D)J `all_contig*`
#' files are renamed with the sample ID as a prefix. Per-sample V(D)J folders
#' are intentionally not used because they contain filtered contig outputs.
#'
#' @param input_dir Character scalar. Project or Cell Ranger output directory to
#'   search.
#' @param dest_dir Character scalar. Destination directory for collected files.
#' @param sample_regex Optional regular expression used to extract sample names.
#'   See [TS_find_cellranger_samples()].
#' @param overwrite Logical. If `FALSE`, existing destination files cause an
#'   error before copying starts.
#' @param confirm Logical. If `TRUE`, ask for confirmation before copying.
#' @param execute Logical. If `FALSE`, return the planned copy operations without
#'   creating directories or copying files.
#' @param strict Logical. If `TRUE`, existing raw V(D)J directories must contain
#'   all requested `all_contig*` files. If `FALSE`, missing V(D)J files are
#'   skipped and recorded in the returned plan.
#' @param ... Additional arguments passed to [TS_find_cellranger_samples()], such
#'   as `matrix_dir_names` or `vdj_dir_names`.
#'
#' @return Invisibly returns the copy plan tibble.
#'
#' @examples
#' project <- file.path(tempdir(), "cellranger_project_collect")
#' matrix_dir <- file.path(project, "sampleA", "cellranger_outs", "filtered_feature_bc_matrix")
#' dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(matrix_dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")))
#'
#' TS_collect_cellranger_files(
#'   input_dir = project,
#'   dest_dir = file.path(tempdir(), "cellranger_transfer_collect"),
#'   execute = FALSE
#' )
#'
#' @export
TS_collect_cellranger_files <- function(
  input_dir,
  dest_dir,
  sample_regex = NULL,
  overwrite = FALSE,
  confirm = interactive(),
  execute = TRUE,
  strict = TRUE,
  ...
) {
  .TS_check_logical_scalar(confirm, "confirm")
  .TS_check_logical_scalar(execute, "execute")

  plan <- TS_plan_cellranger_file_collection(
    input_dir = input_dir,
    dest_dir = dest_dir,
    sample_regex = sample_regex,
    overwrite = overwrite,
    strict = strict,
    ...
  )

  printable_plan <- dplyr::select(
    plan,
    dplyr::all_of(c(
      "sample_id",
      "file_group",
      "source_path",
      "dest_path",
      "status"
    ))
  )

  message("Proposed Cell Ranger file collection:")
  print(printable_plan, n = Inf)

  copy_plan <- dplyr::filter(plan, .data$copy)

  if (nrow(copy_plan) == 0) {
    message("Nothing to copy.")
    return(invisible(plan))
  }

  if (!isTRUE(execute)) {
    message(
      "Dry run only. No directories were created and no files were copied."
    )
    return(invisible(plan))
  }

  if (isTRUE(confirm)) {
    response <- readline(
      prompt = paste0(
        "Proceed with copying ",
        nrow(copy_plan),
        " files? [y/N]: "
      )
    )

    if (!tolower(trimws(response)) %in% c("y", "yes")) {
      message("Aborted. No files were copied.")
      return(invisible(plan))
    }
  }

  unique_dest_dirs <- unique(copy_plan$dest_dir)
  vapply(
    unique_dest_dirs,
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  )

  copy_ok <- file.copy(
    from = copy_plan$source_path,
    to = copy_plan$dest_path,
    overwrite = overwrite,
    copy.date = FALSE
  )

  if (any(!copy_ok)) {
    failed <- copy_plan$dest_path[!copy_ok]
    stop(
      "Some files could not be copied:\n",
      paste(failed, collapse = "\n"),
      call. = FALSE
    )
  }

  message("Files copied successfully.")
  invisible(plan)
}

.TS_check_character_scalar <- function(x, arg) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop(arg, " must be a non-empty character scalar.", call. = FALSE)
  }
}

.TS_check_logical_scalar <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(arg, " must be TRUE or FALSE.", call. = FALSE)
  }
}

.TS_list_dirs <- function(path) {
  dirs <- list.dirs(path, full.names = TRUE, recursive = TRUE)
  normalizePath(dirs, winslash = "/", mustWork = TRUE)
}

.TS_path_parts <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  strsplit(normalized, "/", fixed = TRUE)[[1]]
}

.TS_path_from_parts <- function(parts, index) {
  if (length(parts) > 0 && identical(parts[[1]], "")) {
    paste0("/", paste(parts[seq.int(2, index)], collapse = "/"))
  } else {
    paste(parts[seq_len(index)], collapse = "/")
  }
}

.TS_is_combined_multi_matrix_dir <- function(matrix_dirs) {
  vapply(
    matrix_dirs,
    function(matrix_dir) {
      parts <- .TS_path_parts(matrix_dir)
      has_per_sample_marker <- any(
        parts %in% c("per_sample_outs", "per_samples_outs")
      )

      if (has_per_sample_marker) {
        return(FALSE)
      }

      parent_dir <- dirname(matrix_dir)
      any(dir.exists(file.path(
        parent_dir,
        c("per_sample_outs", "per_samples_outs")
      )))
    },
    logical(1)
  )
}

.TS_extract_regex_sample <- function(path, sample_regex) {
  if (is.null(sample_regex)) {
    return(NA_character_)
  }

  match <- regexec(sample_regex, path, perl = TRUE)
  pieces <- regmatches(path, match)[[1]]

  if (length(pieces) == 0) {
    return(NA_character_)
  }

  if (length(pieces) >= 2) {
    pieces[[2]]
  } else {
    pieces[[1]]
  }
}

.TS_infer_cellranger_sample <- function(matrix_dir, sample_regex = NULL) {
  parts <- .TS_path_parts(matrix_dir)
  per_sample_markers <- which(
    parts %in% c("per_sample_outs", "per_samples_outs")
  )
  regex_sample <- .TS_extract_regex_sample(matrix_dir, sample_regex)
  regex_matched <- !is.na(regex_sample) && nzchar(regex_sample)

  if (length(per_sample_markers) > 0) {
    marker_index <- per_sample_markers[[length(per_sample_markers)]]

    if (length(parts) <= marker_index) {
      stop(
        "Could not infer sample name from per_sample_outs path: ",
        matrix_dir,
        call. = FALSE
      )
    }

    sample_index <- marker_index + 1
    sample_id <- parts[[sample_index]]
    sample_root <- .TS_path_from_parts(parts, sample_index)
    outs_dir <- .TS_path_from_parts(parts, marker_index - 1)
    detection_note <- paste0(parts[[marker_index]], "/<sample>")
  } else {
    output_markers <- which(parts %in% c("cellranger_outs", "outs"))

    if (length(output_markers) > 0 && output_markers[[1]] > 1) {
      output_index <- output_markers[[length(output_markers)]]
      sample_index <- output_index - 1
      sample_id <- parts[[sample_index]]
      sample_root <- .TS_path_from_parts(parts, output_index)
      outs_dir <- sample_root
      detection_note <- paste0(parts[[output_index]], " parent directory")
    } else {
      sample_id <- basename(dirname(matrix_dir))
      sample_root <- dirname(matrix_dir)
      outs_dir <- sample_root
      detection_note <- "matrix parent directory"
    }
  }

  if (regex_matched) {
    sample_id <- regex_sample
    detection_note <- paste0(detection_note, "; sample_regex")
  }

  list(
    sample_id = sample_id,
    sample_root = normalizePath(sample_root, winslash = "/", mustWork = TRUE),
    outs_dir = normalizePath(outs_dir, winslash = "/", mustWork = TRUE),
    detection_note = detection_note
  )
}

.TS_find_raw_vdj_dirs <- function(outs_dir, vdj_dir_names) {
  if (!dir.exists(outs_dir)) {
    return(character())
  }

  direct <- file.path(outs_dir, vdj_dir_names)
  direct <- direct[dir.exists(direct)]
  unique(normalizePath(direct, winslash = "/", mustWork = TRUE))
}

.TS_first_or_na <- function(x) {
  if (length(x) == 0) {
    NA_character_
  } else {
    x[[1]]
  }
}

.TS_validate_sample_ids <- function(sample_ids) {
  bad <- is.na(sample_ids) |
    !nzchar(sample_ids) |
    sample_ids %in% c(".", "..") |
    grepl("[/\\\\]", sample_ids)

  if (any(bad)) {
    stop(
      "Some sample IDs are missing or unsafe for use as folder names:\n",
      paste(sample_ids[bad], collapse = "\n"),
      call. = FALSE
    )
  }
}

.TS_validate_cellranger_vdj_sample_scope <- function(samples) {
  samples_with_vdj <- dplyr::filter(
    samples,
    !is.na(.data$vdj_b_dir) | !is.na(.data$vdj_t_dir)
  )

  if (nrow(samples_with_vdj) == 0) {
    return(invisible(samples))
  }

  outs_counts <- dplyr::summarise(
    dplyr::group_by(samples_with_vdj, .data$outs_dir),
    n = dplyr::n_distinct(.data$sample_id),
    .groups = "drop"
  )
  multi_sample_outs <- dplyr::filter(outs_counts, .data$n > 1)

  if (nrow(multi_sample_outs) == 0) {
    return(invisible(samples))
  }

  affected_samples <- dplyr::filter(
    samples_with_vdj,
    .data$outs_dir %in% multi_sample_outs$outs_dir
  )

  stop(
    "Raw V(D)J all_contig files are run-level outputs, not sample-level outputs. ",
    "Multiple samples were detected under the same outs directory, so these ",
    "raw V(D)J files will not be copied into per-sample folders automatically:\n",
    paste(
      paste0(
        affected_samples$outs_dir,
        " -> ",
        affected_samples$sample_id
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

.TS_plan_one_cellranger_sample <- function(
  sample_row,
  dest_dir,
  matrix_files,
  vdj_files
) {
  sample_id <- sample_row$sample_id[[1]]

  matrix_plan <- .TS_make_file_plan(
    sample_id = sample_id,
    file_group = "filtered_feature_bc_matrix",
    source_dir = sample_row$matrix_dir[[1]],
    dest_dir = file.path(dest_dir, sample_id, "filtered_feature_bc_matrix"),
    source_filenames = matrix_files,
    dest_filenames = matrix_files,
    required = TRUE
  )

  vdj_plans <- lapply(c("vdj_b", "vdj_t"), function(vdj_group) {
    vdj_dir <- sample_row[[paste0(vdj_group, "_dir")]][[1]]

    if (is.na(vdj_dir) || !nzchar(vdj_dir)) {
      return(NULL)
    }

    .TS_make_file_plan(
      sample_id = sample_id,
      file_group = vdj_group,
      source_dir = vdj_dir,
      dest_dir = file.path(dest_dir, sample_id, vdj_group),
      source_filenames = vdj_files,
      dest_filenames = paste0(sample_id, "_", vdj_files),
      required = FALSE
    )
  })

  dplyr::bind_rows(c(list(matrix_plan), vdj_plans))
}

.TS_make_file_plan <- function(
  sample_id,
  file_group,
  source_dir,
  dest_dir,
  source_filenames,
  dest_filenames,
  required
) {
  source_paths <- file.path(source_dir, source_filenames)
  dest_paths <- file.path(dest_dir, dest_filenames)
  source_exists <- file.exists(source_paths)

  tibble::tibble(
    sample_id = sample_id,
    file_group = file_group,
    source_dir = source_dir,
    source_filename = source_filenames,
    source_path = normalizePath(source_paths, winslash = "/", mustWork = FALSE),
    dest_dir = normalizePath(dest_dir, winslash = "/", mustWork = FALSE),
    dest_filename = dest_filenames,
    dest_path = normalizePath(dest_paths, winslash = "/", mustWork = FALSE),
    required = required,
    source_exists = source_exists,
    copy = source_exists,
    status = dplyr::case_when(
      source_exists ~ "ready",
      required ~ "missing_required_source",
      TRUE ~ "missing_optional_source"
    )
  )
}
