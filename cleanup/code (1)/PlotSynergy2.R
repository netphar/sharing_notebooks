PlotSynergy2 = function (data, type = "2D", save.file = FALSE, len = 3, pair.index = NULL, 
          legend.start = NULL, legend.end = NULL, row.range = NULL, 
          col.range = NULL) 
{
  if (!is.list(data)) {
    stop("Input data is not a list format!")
  }
  scores <- data$scores
  drug.pairs <- data$drug.pairs
  num.pairs <- 1:length(scores)
  plots <- list()
  if (!is.null(pair.index)) {
    num.pairs <- pair.index
  }
  for (i in num.pairs) {
    scores.dose <- scores[[i]]
    drug.col <- drug.pairs$drug.col[i]
    drug.row <- drug.pairs$drug.row[i]
    scores.tmp <- scores.dose
    if (!is.null(col.range)) {
      if (col.range[1] == 1) {
        scores.tmp <- scores.tmp[, (col.range[1] + 1):col.range[2]]
      }
      else {
        scores.tmp <- scores.tmp[, col.range[1]:col.range[2]]
      }
    }
    else {
      scores.tmp <- scores.tmp[, -1]
    }
    if (!is.null(row.range)) {
      if (row.range[1] == 1) {
        scores.tmp <- scores.tmp[(row.range[1] + 1):row.range[2], 
                                 ]
      }
      else {
        scores.tmp <- scores.tmp[row.range[1]:row.range[2], 
                                 ]
      }
    }
    else {
      scores.tmp <- scores.tmp[-1, ]
    }
    summary.score <- round(mean(scores.tmp, na.rm = TRUE), 
                           3)
    row.conc <- as.numeric(rownames(scores.dose))
    col.conc <- as.numeric(colnames(scores.dose))
    nr <- nrow(scores.dose)
    nc <- ncol(scores.dose)
    scores.extended <- .ExtendedScores(scores.dose, len)
    mat.tmp <- scores.extended
    if (!is.null(col.range)) {
      mat.tmp[, which(pred.cory >= col.range[1] & pred.cory <= 
                        col.range[2])]
    }
    if (!is.null(row.range)) {
      mat.tmp[which(pred.corx >= row.range[1] & pred.corx <= 
                      row.range[2]), ]
    }
    pp <- matrix(NA, nrow(mat.tmp) * ncol(mat.tmp), 3)
    pp <- data.frame(pp)
    colnames(pp) <- c("x", "y", "z")
    pp$x <- rep(colnames(mat.tmp), each = nrow(mat.tmp))
    pp$y <- rep(rownames(mat.tmp), ncol(mat.tmp))
    pp$z <- c(mat.tmp)
    plot.title <- paste("Average synergy: ", summary.score, 
                        " (", data$method, ")", sep = "")
    conc.runit <- drug.pairs$concRUnit[i]
    conc.cunit <- drug.pairs$concCUnit[i]
    unit.rtext <- paste("(", conc.runit, ")", sep = "")
    unit.ctext <- paste("(", conc.cunit, ")", sep = "")
    file.name <- paste(drug.row, drug.col, "synergy", drug.pairs$blockIDs[i], 
                       data$method, "pdf", sep = ".")
    drug.row <- paste(drug.row, unit.rtext, sep = " ")
    drug.col <- paste(drug.col, unit.ctext, sep = " ")
    max.dose <- max(abs(max(scores.dose)), abs(min(scores.dose)))
    color.range <- round(max.dose + 5, -1)
    if (is.null(legend.start)) {
      start.point <- -color.range
    }
    else {
      start.point <- legend.start
    }
    if (is.null(legend.end)) {
      end.point <- color.range
    }
    else {
      end.point <- legend.end
    }
    levels <- seq(start.point, end.point, by = 2)
    col1 <- colorRampPalette(c("green", "#FFFFFF"))(length(which(levels <= 
                                                                   0)))
    col2 <- colorRampPalette(c("#FFFFFF", "red"))(length(which(levels >= 
                                                                 0)))
    col <- c(col1, col2[-1])
    if (type == "3D") {
      if (!is.null(col.range)) {
        xaxis <- list(at = seq(1, ncol(mat.tmp), by = len + 
                                 1), labels = round(col.conc[col.range[1]:col.range[2]], 
                                                    3))
      }
      else {
        xaxis <- list(at = seq(1, ncol(mat.tmp), by = len + 
                                 1), labels = round(col.conc, 3))
      }
      if (!is.null(row.range)) {
        yaxis <- list(at = seq(1, nrow(mat.tmp), by = len + 
                                 1), labels = round(row.conc[row.range[1]:row.range[2]], 
                                                    3))
      }
      else {
        yaxis <- list(at = seq(1, nrow(mat.tmp), by = len + 
                                 1), labels = round(row.conc, 3))
      }
      par1 <- list(arrows = FALSE, distance = c(0.8, 0.8, 
                                                0.8), col = 1, cex = 0.8, z = list(tick.number = 6), 
                   x = xaxis, y = yaxis)
      zlabs <- list(expression("Synergy score"), rot = 90, 
                    cex = 1, axis.key.padding = 0)
      xpar <- list(as.character(drug.col), cex = 1, rot = 20)
      ypar <- list(as.character(drug.row), cex = 1, rot = -50)
      fig <- wireframe(t(mat.tmp), scales = par1, drape = TRUE, 
                       colorkey = list(space = "top", width = 0.5), 
                       screen = list(z = 30, x = -55), zlab = zlabs, 
                       xlab = xpar, ylab = ypar, zlim = c(start.point, 
                                                          end.point), col.regions = col, main = plot.title, 
                       at = do.breaks(c(start.point, end.point), length(col)), 
                       par.settings = list(axis.line = list(col = "transparent")), 
                       zoom = 1, aspect = 1)
      print(fig)
      fig <- recordPlot()
    }
    else if (type == "2D") {
      dev.new(noRStudioGD = TRUE)
      layout(matrix(c(1, 2), nrow = 2L, ncol = 1L), heights = c(0.1, 
                                                                1))
      par(mar = c(0, 6.1, 2.1, 4.1))
      suppressWarnings(par(mgp = c(3, -1.5, 0)))
      plot.new()
      plot.window(ylim = c(0, 1), xlim = range(levels), 
                  xaxs = "i", yaxs = "i")
      rect(levels[-length(levels)], 0, levels[-1L], 0.3, 
           col = col, border = NA)
      axis(3, tick = FALSE, at = do.breaks(c(start.point, 
                                             end.point), end.point/10))
      title(plot.title)
      par(mar = c(5.1, 4.1, 1.1, 2.1))
      suppressWarnings(par(mgp = c(2, 1, 0)))
      plot.new()
      mat.tmp <- t(mat.tmp)
      x.2D <- (1:dim(mat.tmp)[1] - 1)/(dim(mat.tmp)[1] - 
                                         1)
      y.2D <- (1:dim(mat.tmp)[2] - 1)/(dim(mat.tmp)[2] - 
                                         1)
      plot.window(asp = NA, xlim = range(x.2D), ylim = range(y.2D), 
                  "", xaxs = "i", yaxs = "i")
      .filled.contour(x.2D, y.2D, z = mat.tmp, levels, 
                      col = col)
      box()
      mtext(drug.col, 1, cex = 1, padj = 3)
      mtext(drug.row, 2, cex = 1, las = 3, padj = -3)
      if (!is.null(row.range)) {
        yconc <- round(row.conc[row.range[1]:row.range[2]], 
                       3)
      }
      else {
        yconc <- round(row.conc, 3)
      }
      if (!is.null(col.range)) {
        xconc <- round(col.conc[col.range[1]:col.range[2]], 
                       3)
      }
      else {
        xconc <- round(col.conc, 3)
      }
      axis(side = 1, at = seq(0, 1, by = 1/(length(xconc) - 
                                              1)), labels = xconc)
      axis(side = 2, at = seq(0, 1, by = 1/(length(yconc) - 
                                              1)), labels = yconc)
      fig <- recordPlot()
      dev.off()
    }
    else {
      if (!is.null(col.range)) {
        xaxis <- list(at = seq(1, ncol(mat.tmp), by = len + 
                                 1), labels = round(col.conc[col.range[1]:col.range[2]], 
                                                    3))
      }
      else {
        xaxis <- list(at = seq(1, ncol(mat.tmp), by = len + 
                                 1), labels = round(col.conc, 3))
      }
      if (!is.null(row.range)) {
        yaxis <- list(at = seq(1, nrow(mat.tmp), by = len + 
                                 1), labels = round(row.conc[row.range[1]:row.range[2]], 
                                                    3))
      }
      else {
        yaxis <- list(at = seq(1, nrow(mat.tmp), by = len + 
                                 1), labels = round(row.conc, 3))
      }
      par1 <- list(arrows = FALSE, distance = c(0.8, 0.8, 
                                                0.8), col = 1, cex = 0.8, z = list(tick.number = 6), 
                   x = xaxis, y = yaxis)
      zlabs <- list(expression("Synergy score"), rot = 90, 
                    cex = 1, axis.key.padding = 0)
      xpar <- list(as.character(drug.col), cex = 1, rot = 20)
      ypar <- list(as.character(drug.row), cex = 1, rot = -50)
      syn.3d.plot <- wireframe(t(mat.tmp), scales = par1, 
                               drape = TRUE, colorkey = list(space = "top", 
                                                             width = 0.5), screen = list(z = 30, x = -55), 
                               zlab = zlabs, xlab = xpar, ylab = ypar, zlim = c(start.point, 
                                                                                end.point), col.regions = col, main = plot.title, 
                               at = do.breaks(c(start.point, end.point), length(col)), 
                               par.settings = list(axis.line = list(col = "transparent")), 
                               zoom = 1, aspect = 1)
      layout(matrix(c(1, 2, 3, 3), nrow = 2L, ncol = 2L), 
             heights = c(0.1, 1))
      par(mar = c(0, 6.1, 2.1, 4.1))
      suppressWarnings(par(mgp = c(3, -0.8, 0)))
      plot.new()
      plot.window(ylim = c(0, 1), xlim = range(levels), 
                  xaxs = "i", yaxs = "i")
      rect(levels[-length(levels)], 0, levels[-1L], 0.3, 
           col = col, border = NA)
      axis(3, tick = FALSE, at = do.breaks(c(start.point, 
                                             end.point), end.point/10))
      title(plot.title)
      par(mar = c(5.1, 4.1, 1.1, 2.1))
      suppressWarnings(par(mgp = c(2, 1, 0)))
      plot.new()
      mat.tmp <- t(mat.tmp)
      x.2D <- (1:dim(mat.tmp)[1] - 1)/(dim(mat.tmp)[1] - 
                                         1)
      y.2D <- (1:dim(mat.tmp)[2] - 1)/(dim(mat.tmp)[2] - 
                                         1)
      plot.window(asp = NA, xlim = range(x.2D), ylim = range(y.2D), 
                  "", xaxs = "i", yaxs = "i")
      .filled.contour(x.2D, y.2D, z = mat.tmp, levels, 
                      col = col)
      box()
      mtext(drug.col, 1, cex = 1, padj = 3)
      mtext(drug.row, 2, cex = 1, las = 3, padj = -3)
      if (!is.null(row.range)) {
        yconc <- round(row.conc[row.range[1]:row.range[2]], 
                       3)
      }
      else {
        yconc <- round(row.conc, 3)
      }
      if (!is.null(col.range)) {
        xconc <- round(col.conc[col.range[1]:col.range[2]], 
                       3)
      }
      else {
        xconc <- round(col.conc, 3)
      }
      axis(side = 1, at = seq(0, 1, by = 1/(length(xconc) - 
                                              1)), labels = xconc)
      axis(side = 2, at = seq(0, 1, by = 1/(length(yconc) - 
                                              1)), labels = yconc)
      print(syn.3d.plot, position = c(0.5, 0, 1, 1), newpage = FALSE)
      fig <- recordPlot()
    }
    plots[[i]] <- fig
    if (save.file) {
      if (type %in% c("2D", "3D")) {
        pdf(file.name, width = 10, height = 10)
      }
      else {
        pdf(file.name, width = 12, height = 6)
      }
      print(fig)
      dev.off()
    }
  }
  if (!save.file) {
    for (i in num.pairs) {
      dev.new(noRStudioGD = TRUE)
      replayPlot(plots[[i]])
    }
  }
}
