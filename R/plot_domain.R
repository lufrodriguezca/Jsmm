# R/plot_domain.R

#' Plot domain
#'
#' Create a plot from the polygonal triangulation.
#'
#' @param domain polygonal object from the class sf.
#' @param customize_plot list with optional presetting for visualization.
#' @return A customized domain plot.
#' @export
#'

plot_domain = function(domain, customize_plot = NULL){
  customize_plot = preset_dom_plot(domain, customize_plot)

  center_location = centroid_locations(domain$polygon)
  nsites  = nrow(center_location)

  (ggplot2::ggplot(data = domain$polygon)
    + ggplot2::geom_sf(data = domain$polygon[which(domain$polygon$id > 0), ], aes(fill = factor(id)),
              color = customize_plot$color_polyb, lwd = customize_plot$lwd_poly)
    + scale_fill_manual(values = customize_plot$color_map, labels = customize_plot$lab_legend, name = customize_plot$name_legend)
    + ggplot2::geom_sf(data = domain$polygon[which(domain$polygon$id_location != 0), ], fill = customize_plot$site_col , color = customize_plot$siteb_col,
              lwd = customize_plot$lwd_poly)
    + ggplot2::geom_text(data = center_location, aes(x = x_centroid, y = y_centroid, label = id_location), size = customize_plot$site_text_size,
                nudge_y = rep(customize_plot$site_text_pos[1], nsites),
                nudge_x = rep(customize_plot$site_text_pos[2], nsites))
    + xlab(customize_plot$labs_axis[1])
    + ylab(customize_plot$labs_axis[2])

    + theme(panel.border = element_rect(fill = NA, color = customize_plot$panelb_color, linewidth = customize_plot$panelb_lwd),
            panel.background = element_rect(fill = NA))
  )
}
