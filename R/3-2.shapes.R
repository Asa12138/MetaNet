# ========3.1.shapes========

add_metanet_shape_diamond <- function() {
  mydiamond <- function(coords, v = NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]
    }
    vertex.size <- 1 / 200 * sqrt(2) * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    vertex.size <- rep(vertex.size, length.out = nrow(coords))
    vertex.frame.color[vertex.frame.width <= 0] <- NA
    vertex.frame.width[vertex.frame.width <= 0] <- 1
    if (length(vertex.frame.width) == 1) {
      symbols(
        x = coords[, 1], y = coords[, 2], bg = vertex.color,
        fg = vertex.frame.color, stars = cbind(vertex.size, vertex.size, vertex.size, vertex.size),
        lwd = vertex.frame.width, add = TRUE, inches = FALSE
      )
    } else {
      mapply(coords[, 1], coords[, 2], vertex.color, vertex.frame.color,
        vertex.size, vertex.frame.width,
        FUN = function(x, y, bg, fg, size, lwd) {
          symbols(
            x = x, y = y, bg = bg, fg = fg, lwd = lwd,
            stars = cbind(size, size, size, size), add = TRUE, inches = FALSE
          )
        }
      )
    }
  }
  igraph::add_shape("diamond", clip = shape_noclip, plot = mydiamond)
}

add_metanet_shape_triangle1 <- function() {
  mytriangle1 <- function(coords, v = NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]
    }
    vertex.size <- 1 / 200 * 1.2 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    vertex.size <- rep(vertex.size, length.out = nrow(coords))
    vertex.frame.color[vertex.frame.width <= 0] <- NA
    vertex.frame.width[vertex.frame.width <= 0] <- 1
    if (length(vertex.frame.width) == 1) {
      symbols(
        x = coords[, 1], y = coords[, 2] - vertex.size / sqrt(3), bg = vertex.color,
        fg = vertex.frame.color, stars = cbind(vertex.size, vertex.size * sqrt(3), vertex.size, 0),
        lwd = vertex.frame.width, add = TRUE, inches = FALSE
      )
    } else {
      mapply(coords[, 1], coords[, 2], vertex.color, vertex.frame.color,
        vertex.size, vertex.frame.width,
        FUN = function(x, y, bg, fg, size, lwd) {
          symbols(
            x = x, y = y - size / sqrt(3), bg = bg, fg = fg, lwd = lwd,
            stars = cbind(size, size * sqrt(3), size, 0), add = TRUE, inches = FALSE
          )
        }
      )
    }
  }
  igraph::add_shape("triangle1", clip = shape_noclip, plot = mytriangle1)
}

add_metanet_shape_triangle2 <- function() {
  mytriangle2 <- function(coords, v = NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]
    }
    vertex.size <- 1 / 200 * 1.2 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    vertex.size <- rep(vertex.size, length.out = nrow(coords))
    vertex.frame.color[vertex.frame.width <= 0] <- NA
    vertex.frame.width[vertex.frame.width <= 0] <- 1
    if (length(vertex.frame.width) == 1) {
      symbols(
        x = coords[, 1], y = coords[, 2] + vertex.size / sqrt(3), bg = vertex.color,
        fg = vertex.frame.color, stars = cbind(vertex.size, 0, vertex.size, vertex.size * sqrt(3)),
        lwd = vertex.frame.width, add = TRUE, inches = FALSE
      )
    } else {
      mapply(coords[, 1], coords[, 2], vertex.color, vertex.frame.color,
        vertex.size, vertex.frame.width,
        FUN = function(x, y, bg, fg, size, lwd) {
          symbols(
            x = x, y = y + size / sqrt(3), bg = bg, fg = fg, lwd = lwd,
            stars = cbind(size, 0, size, size * sqrt(3)), add = TRUE, inches = FALSE
          )
        }
      )
    }
  }
  igraph::add_shape("triangle2", clip = shape_noclip, plot = mytriangle2)
}

# 因为igraph默认只有circle和square适合展示，所以这里添加了更多的形状

add_metanet_shapes <- function() {
  # ！！！考虑添加更多的形状，至少是21:25，形状为pie时添加额外legend
  for (i in names(default_v_shape)[3:5]) {
    paste0("add_metanet_shape_", i) -> fun_name
    do.call(fun_name, list())
  }
}
