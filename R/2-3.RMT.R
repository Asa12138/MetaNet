# =========2.2RMT optimize=====

#' Get RMT threshold for a correlation matrix
#'
#' @param occor.r a corr object or a correlation matrix
#' @param min_threshold min_threshold
#' @param max_threshold max_threshold
#' @param step step
#' @param gif render a .gif file?
#' @param verbose verbose
#' @param out_dir output dir
#'
#' @return a r-threshold
#' @export
#' @references
#' J. Zhou, Y. Deng, FALSE. Luo, Z. He, Q. Tu, X. Zhi, (2010) Functional Molecular Ecological Networks, doi:10.1128/mBio.00169-10.
#' <https://matstat.org/content_en/RMT/RMThreshold_Intro.pdf>
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' t(otutab) -> totu
#' c_net_calculate(totu) -> corr
#' rmt(corr)
#' # recommend: 0.69
#' c_net_build(corr, r_threshold = 0.69) -> co_net_rmt
#' }
RMT_threshold <- function(occor.r, out_dir, min_threshold = 0.5, max_threshold = 0.8,
                          step = 0.02, gif = FALSE, verbose = FALSE) {
  nwd <- getwd()
  on.exit(setwd(nwd))

  setwd(out_dir)
  if (inherits(occor.r, "corr")) occor.r <- occor.r$r
  if (!dir.exists("./RMT_temp")) dir.create("./RMT_temp")
  diag(occor.r) <- 0

  if (max_threshold >= max(abs(occor.r))) max_threshold <- (max(abs(occor.r)) - step)
  if (min_threshold >= max_threshold) min_threshold <- max_threshold - 10 * step

  thres_seq <- seq(min_threshold, max_threshold, step)

  res <- data.frame()
  for (i in seq_len(length(thres_seq))) {
    threshold <- thres_seq[i]
    if (!verbose) pcutils::dabiao(paste0("Calculating", i, ":  threshold =", signif(threshold, 3)), print = TRUE)
    corr_r1 <- occor.r
    corr_r1[abs(corr_r1) < threshold] <- 0
    # calculate eigenvalues
    rand.mat <- corr_r1
    eigenvalues <- eigen(rand.mat, only.values = TRUE)$values
    eigenvalues <- eigenvalues[order(eigenvalues)] / max(abs(eigenvalues))
    eigenvalues <- pcutils::remove.outliers(unique(eigenvalues))

    # get the NNDS
    { # uf <- rm.unfold.gauss(eigenvalues,pop.up = TRUE)
      dens <- density(eigenvalues, kernel = "gaussian")
      midpoints <- \(x)(x[-length(x)] + 0.5 * diff(x))
      scale.function <- approx(dens$x, dens$y, xout = midpoints(eigenvalues))
      ev.spacing <- diff(eigenvalues)
      ev.spacing <- ev.spacing * scale.function$y
      ev.spacing <- ev.spacing / mean(ev.spacing)
    }

    ev.spacing <- ev.spacing[ev.spacing <= 3]
    # test whether fit possion?
    p_ks_test <- ks.test(unique(ev.spacing), "pexp", 1)$p.value
    # get sse
    # sse = rm.sse(ev.spacing)
    sse <- get_sse(ev.spacing)
    log_sse <- log(sse)

    # maximum likelihood
    evs <- ev.spacing[ev.spacing != 0]
    N <- length(evs)
    log_LE <- -sum(evs) / N
    log_LW <- log(pi / 2) + sum(log(evs)) / N - 0.25 * pi * sum(evs^2) / N

    # save png
    {
      histo <- hist(ev.spacing, breaks = seq(min(ev.spacing), max(ev.spacing), len = 51), plot = FALSE)
      grDevices::png(paste0("RMT_temp/rmt_nnsd", i, ".png"), height = 600, width = 700, res = 130)
      nnsd_plot(
        histo = histo, title = "Eigenvalue spacing distribution (NNSD)", threshold = threshold,
        dis_GOE = log_LW, dis_possion = log_LE, p_ks_test = p_ks_test
      )
      grDevices::dev.off()
    }
    res <- rbind(res, data.frame(threshold, p_ks_test, log_sse, log_LW, log_LE))
  }
  message("The Intermediate files saved in ", out_dir, "/RMT_temp/.")
  # transfer to gif
  if (gif) {
    lib_ps("gifski", library = FALSE)
    gifski::gifski(paste0("RMT_temp/rmt_nnsd", seq_len(length(thres_seq)), ".png"),
      gif_file = "RMT_temp/rmt_nnsd.gif"
    )
  }
  r_threshold <- (res[which(res$log_LW == min(res$log_LW)), "threshold"] +
    res[which(res$log_LE == max(res$log_LE)), "threshold"]) / 2
  res <- list(res = res, r_threshold = r_threshold)
  class(res) <- c("rmt_res", class(res))
  return(res)
}

#' Plot a rmt_res
#'
#' @param x rmt_res
#' @param ... Additional arguments
#'
#' @return ggplot
#' @exportS3Method
#' @method plot rmt_res
plot.rmt_res <- function(x, ...) {
  threshold <- value <- variable <- xi <- y <- NULL
  res <- x$res
  linedf <- data.frame(
    variable = c("p_ks_test", "log_sse", "log_LW", "log_LE"),
    xi = c(
      res[which(res$p_ks_test == max(res$p_ks_test))[1], "threshold"],
      res[which(res$log_sse == min(res$log_sse))[1], "threshold"],
      res[which(res$log_LW == min(res$log_LW))[1], "threshold"],
      res[which(res$log_LE == max(res$log_LE))[1], "threshold"]
    ),
    x = max(res$threshold) - min(res$threshold),
    y = apply(res[, -1], 2, max)
  )

  reshape2::melt(res, "threshold") -> md

  # filter(threshold<0.77)%>%
  p <- ggplot(md, aes(threshold, value)) +
    geom_point(aes(col = variable)) +
    geom_line(aes(col = variable)) +
    scale_color_manual(values = get_cols(4, "col1")) +
    facet_wrap(. ~ variable, scales = "free_y") +
    theme_bw() +
    xlab(NULL) +
    geom_text(data = linedf, aes(x = xi - 0.02 * x, y = 0.5 * y, label = xi)) +
    geom_vline(data = linedf, aes(xintercept = xi), linetype = 2, col = "red") +
    theme(legend.position = "none")

  message(paste("recommend r_threshold: ", mean(linedf$xi)))
  return(p)
}

nnsd_plot <- \(histo = histo, title = title, threshold = threshold,
  dis_GOE = dis_GOE, dis_possion = dis_possion, p_ks_test = p_ks_test) {
  plot(histo, freq = FALSE, col = "#F4FCA1", main = title, font.main = 1, xlab = "eigenvalue spacing", ylab = "PDF of eigenvalue spacing")
  {
    actual.ymax <- par("yaxp")[2]
    x0 <- -log(actual.ymax * 0.98)
    possion_dis <- \(x)exp(-x)
    graphics::curve(possion_dis,
      from = max(x0, min(histo$breaks)),
      to = max(histo$breaks), n = 1001, add = TRUE, type = "l", lty = 1, col = "#EB34FF", lwd = 2
    )
  }
  {
    GOE <- function(x) pi / 2 * x * exp(-pi / 4 * x^2)
    graphics::curve(GOE,
      from = min(histo$breaks),
      to = max(histo$breaks), n = 1001, add = TRUE, type = "l",
      lty = 1, col = "blue", lwd = 2
    )
  }

  if ((!is.na(dis_GOE)) && (!is.na(dis_possion))) {
    graphics::mtext(side = 3, paste(
      "Distance to GOE =", signif(dis_GOE, 3),
      "\nDistance to Possion =", signif(dis_possion, 3), "; ks_test p.value for possion =", signif(p_ks_test, 3)
    ), col = "#878787", cex = 0.6)
  }

  if (!is.na(threshold)) graphics::mtext(side = 4, paste("threshold =", signif(threshold, 4)))

  graphics::legend("topright", inset = 0.05, c("Possion", "GOE"), col = c("#EB34FF", "blue"), lty = 1, lwd = 2, cex = 0.8)
}

trapez <- \(x, y){
  ind <- 2:length(x)
  as.double((x[ind] - x[ind - 1]) %*% (y[ind] + y[ind - 1])) / 2
}

get_sse <- \(ev.spacing){
  dens <- density(ev.spacing)
  N <- 20
  x <- seq(min(ev.spacing), max(ev.spacing), len = 1000)
  A <- exp(-min(ev.spacing)) - exp(-max(ev.spacing))
  xs <- numeric(N + 1)
  xs[1] <- min(ev.spacing)
  for (i in 1:N) xs[i + 1] <- -log(exp(-xs[i]) - A / N)
  area <- numeric(N)
  for (i in 1:N) {
    xsec <- x[(x > xs[i]) & (x < xs[i + 1])]
    xsec <- c(xs[i], xsec, xs[i + 1])
    ysec <- approx(dens$x, dens$y, xout = xsec)$y
    area[i] <- trapez(xsec, ysec)
  }
  sse <- sum((area[i] - A / N)^2)
  sse
}

#' Get RMT threshold for a correlation matrix roughly
#'
#' @export
#' @return recommend threshold
#' @rdname RMT_threshold
rmt <- function(occor.r, min_threshold = 0.5, max_threshold = 0.85, step = 0.01) {
  if (inherits(occor.r, "corr")) occor.r <- occor.r$r
  NNSD <- \(x)abs(diff(x))

  s <- seq(0, 3, 0.1)
  poisson_d <- exp(-s)
  nnsdpois <- density(NNSD(poisson_d))

  ps <- c()
  threshold <- c()

  for (i in seq(min_threshold, max_threshold, step)) {
    corr_r1 <- occor.r
    corr_r1[abs(corr_r1) < i] <- 0
    {
      eigen_res <- sort(eigen(corr_r1)$value)
      # spline to eigen_res
      check <- tryCatch(ssp <- smooth.spline(eigen_res, control.spar = list(low = 0, high = 3)),
        error = \(e) {
          TRUE
        }
      )
      if (rlang::is_true(check)) next
      nnsdw <- density(NNSD(ssp$y))
      chival <- sum((nnsdw$y - nnsdpois$y)^2 / nnsdpois$y / 1e3)
    }

    ps <- c(ps, chival)
    threshold <- c(threshold, i)
    if (((i * 100) %% 5 == 0)) {
      message(paste0("Calculating: ", i))
    }
  }

  res <- data.frame(threshold, ps)
  recommend_thres <- res[which.min(res[, 2]), 1]
  p <- ggplot(res, aes(threshold, ps)) +
    geom_point() +
    geom_vline(xintercept = recommend_thres, linetype = 2, col = "red") +
    geom_text(x = recommend_thres + 0.01, y = 0.5 * max(res$ps), label = recommend_thres) +
    theme_bw(base_size = 15)
  print(p)

  res1 <- res[(res$threshold < (recommend_thres + 0.05)) & (res$threshold > (recommend_thres - 0.05)), ]
  p <- ggplot(res1, aes(threshold, ps)) +
    geom_point() +
    geom_vline(xintercept = recommend_thres, linetype = 2, col = "red") +
    geom_text(x = recommend_thres + 0.01, y = 0.5 * max(res1$ps), label = recommend_thres) +
    theme_bw(base_size = 15)
  print(p)

  message("We recommend r-threshold: ", recommend_thres, ", you can calculate again in a smaller region")
  recommend_thres
}
