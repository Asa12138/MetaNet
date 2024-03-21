#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Assignment pipe
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%<>\%}} for details.
#'
#' @name %<>%
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

# flag determined by the correlation table from one table or two tables
t_flag <- \(corr){
  if (!nrow(corr) == ncol(corr)) {
    return(FALSE)
  }
  if (!all(rownames(corr) == colnames(corr))) {
    return(FALSE)
  }
  return(TRUE)
}

# 多列表是否一一对应
# test 对应的列需要跟剩下的列具有唯一映射
e_match <- function(df, test = NULL) {
  dis_df <- dplyr::distinct_all(df) %>% as.data.frame()
  if (is.null(test)) {
    return(all_same(c(nrow(dis_df), apply(dis_df, 2, \(i)length(unique(i))))))
  } else {
    return(all_same(c(nrow(dis_df), nrow(dplyr::distinct_all(dis_df[, test, drop = FALSE])))))
  }
}

all_same <- \(x){
  return(all(x == x[1]))
}

# choose the last not na value
condance <- \(aa){
  aa <- as.data.frame(aa)
  if (any(is.na(aa[, length(aa)]))) {
    res <- apply(aa, 1, \(x){
      tmp <- x[!is.na(x)]
      if (length(tmp) == 0) {
        return(NA)
      }
      return(tmp[length(tmp)])
    })
  } else {
    res <- aa[, length(aa), drop = TRUE]
  }
  res
}

# bind two df with same columns, the last df will replace first df
cbind_new <- \(df, df1){
  if (ncol(df) < 1) {
    return(df1)
  }
  if (ncol(df1) < 1) {
    return(df)
  }
  inter <- intersect(colnames(df1), colnames(df))
  la <- setdiff(colnames(df), inter)
  cbind(df[, la, drop = FALSE], df1)
}

deprecated <- function(old, new) {
  assign(old, new, envir = asNamespace(packageName()))
}

#' @export c_net_cal
deprecated("c_net_cal", c_net_calculate)

#' @export c_net_module
deprecated("c_net_module", module_detect)

#' @export c_net_lay
deprecated("c_net_lay", c_net_layout)

#' @export as.metanet
deprecated("as.metanet", c_net_update)
