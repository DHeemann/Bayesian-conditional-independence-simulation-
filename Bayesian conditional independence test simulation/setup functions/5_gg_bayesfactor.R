bf_segments <- data.frame(xmin = -Inf,
                          xmax = Inf,
                          ymin = c(1/rev(c(3, 10, 30, 100)),
                                   c(1, 3, 10, 30, 100))) %>%
  mutate(ymax = lead(ymin)) %>%
  rowwise() %>%
  mutate(ymax = replace_na(ymax, ymin*3)) %>%
  ungroup() %>%
  mutate(col = c("grey50","grey70","grey80",
                 "grey100","grey100","grey80",
                 "grey70","grey50","grey30")) 
                          
                                   

theme_bayesfactor <- ggplot() +
  geom_rect(data = bf_segments, 
            aes(xmin = -Inf, 
                xmax = Inf, 
                ymin = ymin, 
                ymax = ymax,
                fill =  col),
            alpha = 0.3) +
  scale_y_continuous(trans='log10', 
                     breaks = c(1/c(3, 10, 30, 100),
                                c(1, 3, 10, 30, 100)),
                     labels = c(paste0(1, "/", c(3, 10, 30, 100)), 
                                paste0(c(1, 3, 10, 30, 100))),
                     sec.axis = sec_axis(trans=~., name="",
                                         breaks = c(1/18, 1/5.5, 1, 5.5,18),
                                         labels = c("strong\nevidence",
                                                    "moderate\nevidence",
                                                    "anectodal\nevidence",
                                                    "moderate\nevidence",
                                                    "strong\nevidence"))) +
  coord_cartesian(ylim = c(1/30, 30))+
  theme_pubr() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_discrete(name = unname(TeX("u_{12}"))) +
  scale_fill_manual(values = c("grey50" = "grey50",
                               "grey70" ="grey70",
                               "grey80"= "grey85",
                               "grey100" = "grey100",
                               "grey30" = "grey30"),
                    guide = FALSE) +
  theme(axis.ticks.y.right = element_blank(),
        #axis.text.x = element_text(angle = -45),
        axis.line.y.right = element_blank(),
        axis.text = element_text(size = 8))


