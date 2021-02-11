blur_point <- function(x, levels = seq(1e-3, 1-1e-1, l = 15),
                       alpha_mult = 1, col = "black", ...){
  for(level in levels){
    polygon(ellipse::ellipse(x, level = level, ...), border = NA, 
            col = scales::alpha(col, min(1, 1/length(levels)*alpha_mult/sqrt(prod(diag(x))))))
  }
}

plot(0, 0, xlim = c(-2, 2), ylim = c(-2, 2), type = "n")
blur_point(0.2 * diag(2), centre = c(1.2, 1.3), alpha_mult = 0.1)
