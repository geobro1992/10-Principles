library(ggplot2)
library(ggrepel)


########################
# amphibians

dat = read.csv("vtparity.csv")
amps = dat[which(dat$class == "Amphibia"),]


p1 <- ggplot(amps, aes(x= juvMort, y = adMort)) + 
  xlim(0, 12.5) +
  ylim(0, 0.9) +
  xlab("Juvenile Mortality") + ylab("Adult Mortality") +
  geom_point(aes(shape = order), size = 2, alpha = 1,
             color = dplyr::case_when(amps$species == "Lithobates catesbeiana" ~ "#1b9e77", 
                                      amps$species == "Plethodon jordanii" ~ "#d95f02",
                                      amps$species == "Anomaloglossus verbeeksnyderorum" ~ "#7570b3",
                                      TRUE ~ "dark grey"))

p1 = p1 + geom_text_repel(aes(label = species), data = subset(amps, species == "Lithobates catesbeiana"),
                  size          = 4,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x",
                  nudge_y       = -0.1,
                  nudge_x       = 1) +
  geom_text_repel(aes(label = species), data = subset(amps, species == "Plethodon jordanii"),
                  size          = 4,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x",
                  nudge_y       = 0,
                  nudge_x       = -3) +
  geom_text_repel(aes(label = species), data = subset(amps, species == "Anomaloglossus verbeeksnyderorum"),
                  size          = 4,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x",
                  nudge_y       = 0.15,
                  nudge_x       = -1) +
  theme_Publication()

p1

ggsave(p1, filename = "amps_scatter.tiff", width = 8, height = 6, dpi = 600)




########################
# amphibians

dat = read.csv("vtparity.csv")
reps = dat[which(dat$class == "Reptilia"),]


p1 <- ggplot(reps, aes(x= juvMort, y = adMort)) + 
  xlim(0, 12.5) + ylim(0, 0.9) +
  xlab("Juvenile Mortality") + ylab("Adult Mortality") +
  geom_point(aes(shape = order), size = 2, alpha = 1,
             color = dplyr::case_when(reps$species == "Vipera aspis" ~ "#1b9e77", 
                                      reps$species == "Dermochelys coriacea" ~ "#d95f02",
                                      reps$species == "Chamaeleo chamaeleon" ~ "#7570b3",
                                      TRUE ~ "dark grey"))

p1 = p1 + geom_text_repel(aes(label = species), data = subset(reps, species == "Vipera aspis"),
                          size          = 4,
                          box.padding   = 1.5,
                          point.padding = 0.5,
                          segment.size  = 0.2,
                          segment.color = "grey50",
                          direction     = "x",
                          nudge_y       = -0.1,
                          nudge_x       = -0.1) +
  geom_text_repel(aes(label = species), data = subset(reps, species == "Dermochelys coriacea"),
                  size          = 4,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x",
                  nudge_y       = 0.1,
                  nudge_x       = 1) +
  geom_text_repel(aes(label = species), data = subset(reps, species == "Chamaeleo chamaeleon"),
                  size          = 4,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x",
                  nudge_y       = 0.4,
                  nudge_x       = 1) +
  theme_Publication()

p1

ggsave(p1, filename = "reps_scatter.tiff", width = 8, height = 6, dpi = 600)
