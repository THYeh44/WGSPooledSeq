#' Generate a figure of  total mutations (left) and normalized mutation density (right)
#'
#' This function generate a figure including total mutations (left) and normalized mutation density (right).
#' @param 1) Input the .txt file generated from pooled_analysis function, 2) Sepcifiy the chromosome, 3) Left max of the axis, 4) Left min of the axis, 5) Left threthold to display gene names, 6) Right min, 7) Right min, and 8) Right threshold.
#' @keywords Mutation density fugure
#' @export
#' @examples
#' mut_fig()

mut_fig <- function(input, chr, left_max, left_min, left_threshold, right_max, right_min, right_threshold){
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  #Left part
  Left <-input %>% filter(NonSynMut >= 1) %>% filter(CHR == chr) %>%
    mutate(NonSynMut_log10 = log10(NonSynMut)) %>%
    ggplot(., aes(x=NonSynMut_log10, y=POS, col=NonSynMut_log10 >= left_threshold)) +
    scale_color_brewer(palette = "Dark2")+
    geom_point(size = 4, alpha = 0.75) +
    coord_cartesian(clip = "off")+
    #scale_x_continuous(limits=c(0,3))+
    scale_x_reverse(limits=c(left_max, 0))+
    scale_y_continuous(position = "right")+
    ggtitle(expression("log"[10]*"Nonsynonymous"))+
    theme_minimal()+
    theme(plot.title = element_text(size=15,hjust=0.5),
          text = element_text(size=15),
          legend.position="none",
          aspect.ratio=8/1,
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          #plot.margin = unit(c(0,0,0,0), "mm"),
          plot.margin = margin(t=0,r=0,b=0,l=0, unit="mm")
          #axis.title.x = element_text(size=8, hjust=0.5)
    )
  #Right part
  Right <- input %>% filter(NonSynMut >= 1) %>% filter(CHR == chr) %>%
    ggplot(., aes(x=NonSynMut_perKb, y=POS, col=NonSynMut_perKb >= right_threshold)) +
    scale_color_brewer(palette = "Dark2")+
    geom_point(size = 4, alpha = 0.75) +
    coord_cartesian(clip = "off")+
    scale_x_continuous(limits=c(0,right_max))+
    ggtitle("Nonsynonymous/kb")+
    #scale_x_reverse(limits=c(3, 0))+
    #coord_flip()+
    theme_minimal()+
    theme(plot.title = element_text(size=15,hjust=0.5),
          text = element_text(size=15),
          legend.position="none",
          aspect.ratio=8/1,
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          #plot.margin = unit(c(0,0,0,-50), "mm"),
          plot.margin = margin(t=0,r=0,b=0,l=-0, unit="mm")
          #axis.title.x = element_text(size=8, hjust=0.99)
    )
  patchwork<-Left + Right + plot_annotation(theme= theme(plot.margin = margin(unit(c(0,0,0,0), "mm"))))
  #plot_layout(guides='collect')
  patchwork +
    plot_annotation(title = chr, theme= theme(plot.title = element_text(hjust=0.5),
                                              plot.margin = margin(unit(c(0,0,0,0), "mm"))))
}
