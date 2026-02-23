# R/preset_dom_plot.R

#' @noRd
#'

preset_dom_plot = function(domain, customize_plot = NULL){
  customize_plot_aux = list()
  nsub = length(unique(domain$polygon$id))
  colors_def = diverge_hcl(100, h = c(180, 70), c = 70, l = c(90, 95))
  customize_plot_aux$color_map = colors_def[50:(50 + (nsub - 1))]
  customize_plot_aux$labs_axis = c("X", "Y")
  customize_plot_aux$lab_legend = NULL
  customize_plot_aux$name_legend = NULL
  customize_plot_aux$color_polyb = "black"
  customize_plot_aux$lwd_poly = 0.01
  customize_plot_aux$site_col = "red"
  customize_plot_aux$siteb_col = "darkred"
  customize_plot_aux$site_text_pos = c(0, 0)
  customize_plot_aux$site_text_size = 5
  customize_plot_aux$panelb_color = "black"
  customize_plot_aux$panelb_lwd = 1

  if(is.null(customize_plot)){
    customize_plot = customize_plot_aux
  }else{
    if(is.null(customize_plot$color_map)){
      customize_plot$color_map =  customize_plot_aux$color_map
    }
    if(is.null(customize_plot$color_polyb)){
      customize_plot$color_polyb = customize_plot_aux$color_polyb
    }
    if(is.null(customize_plot$lwd_poly)){
      customize_plot$lwd_poly = customize_plot_aux$lwd_poly
    }
    if(is.null(customize_plot$site_col)){
      customize_plot$site_col = customize_plot_aux$site_col
    }
    if(is.null(customize_plot$siteb_col)){
      customize_plot$siteb_col = customize_plot_aux$siteb_col
    }
    if(is.null(customize_plot$site_text_pos)){
      customize_plot$site_text_pos =  customize_plot_aux$site_text_pos
    }
    if(is.null(customize_plot$site_text_size)){
      customize_plot$site_text_size = customize_plot_aux$site_text_size
    }
    if(is.null(customize_plot$panelb_color)){
      customize_plot$panelb_color =  customize_plot_aux$panelb_color
    }
    if(is.null(customize_plot$panelb_lwd)){
      customize_plot$panelb_lwd = customize_plot_aux$panelb_lwd
    }
  }

  return(customize_plot)
}
