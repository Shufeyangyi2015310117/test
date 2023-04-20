#' theme_Publication.
#'
#' @description
#' theme_Publication
#'
#' @details theme_Publication
#' @param base_size the base size of text.
#' @param base_family the base type of plot.
#' @return a ggplot.
#' @export
theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line = element_line(colour="black"),
            axis.ticks.x= element_blank(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.text = element_text(size=28),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.spacing  = unit(0, "cm"),
            legend.title = element_blank(),
            plot.margin=unit(c(10,5,5,0),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}
