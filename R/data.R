#' Human IGHV genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGHV genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 159
#' @usage as.factor(string, levels = TS::IGHV_levels)
"IGHV_levels"


#' Human IGHV families in ascending numeric order
#'
#' A vector that can be used to select and arrange IGHV families.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 7
#' @usage as.factor(string, levels = TS::IGHV_family_levels)
"IGHV_family_levels"


#' Human IGHV families in ascending numeric order with IGHV6 and 7 grouped
#'
#' A vector that can be used to select and arrange IGHV families while IGHV6 and 7 are grouped as they only contain one gene each.
#' This requires that IGHV_family is set to IGHV6/7 for those families, which is done automatically by TS_plot_IG_usage(compare_what = "IGHV_family_grouped")
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 6
#' @usage as.factor(string, levels = TS::IGHV_family_levels_grouped)
"IGHV_family_levels_grouped"


#' Human IGHD genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGHV genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 27
#' @usage as.factor(string, levels = TS::IGHD_levels)
"IGHD_levels"


#' Human IGHD families in ascending numeric order
#'
#' A vector that can be used to select and arrange IGHD families.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 7
#' @usage as.factor(string, levels = TS::IGHD_family_levels)
"IGHD_family_levels"


#' Human IGHJ genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGHJ genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 6
#' @usage as.factor(string, levels = TS::IGHJ_levels)
"IGHJ_levels"


#' Human IGKV genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGKV genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 76
#' @usage as.factor(string, levels = TS::IGKV_levels)
"IGKV_levels"


#' Human IGKV families in ascending numeric order
#'
#' A vector that can be used to select and arrange IGKV families.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 7
#' @usage as.factor(string, levels = TS::IGKV_family_levels)
"IGKV_family_levels"


#' Human IGKJ genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGKJ genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 5
#' @usage as.factor(string, levels = TS::IGKJ_levels)
"IGKJ_levels"


#' Human IGLV genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGLV genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 74
#' @usage as.factor(string, levels = TS::IGLV_levels)
"IGLV_levels"


#' Human IGLV families in ascending numeric order
#'
#' A vector that can be used to select and arrange IGLV families.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 11
#' @usage as.factor(string, levels = TS::IGLV_family_levels)
"IGLV_family_levels"


#' Human IGLJ genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGLJ genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 7
#' @usage as.factor(string, levels = TS::IGLJ_levels)
"IGLJ_levels"


#' Human IGLC genes according to their relative locus location
#'
#' A vector that can be used to select and arrange IGLC genes in the correct order.
#' Data obtained from IMGT
#'
#' @docType data
#' @format A character vector of length 7
#' @usage as.factor(string, levels = TS::IGLC_levels)
"IGLC_levels"
