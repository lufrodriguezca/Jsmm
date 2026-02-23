# R/plot_observation_effort.R

#' Plot Observation effort
#'
#' Visualize where and when captures and releases occurred.
#'
#' @param observation_effort Named list with releases and captures information.
#' @param xlab Customized labels for the x-axis.
#' @param ylab Customized labels for the y-axis.
#' @param by_x Ticks in x frequency.
#' @param cex.axisx Size text x.
#' @param cex.axisy Size text y.
#' @param lasy Rotation label y.
#' @param tcky Control tick marks.
#'
#' @return A plot with releases and captures information.
#' @export
#'
plot_observation_effort = function(observation_effort, xlab, ylab, by_x, cex.axisx, cex.axisy, lasy, tcky){

  locationsp = sort(unique(c(observation_effort[["captures"]][["location"]], observation_effort[["releases"]][["location"]])))
  timesp = sort(unique(c(observation_effort[["releases"]][["time"]],  unlist(observation_effort[["captures"]][["time"]] ))))

  xlim = c(min(timesp) - 0.3, max(timesp) + 0.3)
  ylim = c(0, length(observation_effort[["locations"]]))

  if(observation_effort[["method"]]=="CCP"){
    plot(observation_effort[["captures"]][["time"]][1, 1], observation_effort[["captures"]][["time"]][1, 2],
        col = "white", xlab = "", ylab = "", xlim = xlim, ylim = ylim,
        xaxt = "n", yaxt = "n", frame = FALSE)

    for(i in 1:length(observation_effort[["captures"]][["location"]])){
      segments(x0 = observation_effort[["captures"]][["time"]][i, 1], y0 = observation_effort[["captures"]][["location"]][i],
              x1 = observation_effort[["captures"]][["time"]][i, 2], y1 = observation_effort[["captures"]][["location"]][i],
              col="red", lwd = 3)
      points(x = observation_effort[["captures"]][["time"]][i, 1],
            y = observation_effort[["captures"]][["location"]][i],
            col = "red", pch = 16, cex = 0.8)
      points(x = observation_effort[["captures"]][["time"]][i, 2],
            y = observation_effort[["captures"]][["location"]][i],
            bg = "red",col = "black", lwd = 0.9, pch = 21, cex = 0.8)
    }
  }

  if(observation_effort[["method"]] == "ICP"){
    plot(observation_effort[["captures"]][["time"]][1], observation_effort[["captures"]][["time"]][2],
    col = "white", xlab = "", ylab = "", xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n", frame = FALSE)

    for(i in 1:length(observation_effort[["captures"]][["location"]])){
      points(x = observation_effort[["captures"]][["time"]][i],
              y = observation_effort[["captures"]][["location"]][i],
              col = "red", pch =16, cex = 0.8)
    }
  }

  for(i in 1:length(observation_effort[["releases"]][["location"]])){
    points(x = observation_effort[["releases"]][["time"]][i],
          y = observation_effort[["releases"]][["location"]][i] - 0.3,
          bg = "steelblue1", pch = 21, cex = 0.8, col = "black", lwd = 0.9)
  }

  title( xlab = xlab, ylab = ylab)

  axis(side = 2, at = locationsp , labels = locationsp , las = lasy, tck = tcky, lwd = 0, cex.axis = cex.axisy)
  axis(side = 1, at = seq(min(timesp), max(timesp), by = by_x) , labels = seq(min(timesp), max(timesp), by = by_x), cex.axis = cex.axisx)
}
