install.packages("tidyverse")
library(tidyverse)
install.packages("ggsignif")
library(ggsignif)
install.packages("ggpubr")
library(ggpubr)

#Set up Normalized Control and Sleep Data
ctrl <- c(0.330649011, 0.346783364, 0.328704873)

sleep <- c(0.380951092, 0.366219174, 0.402398523, 0.407813692, 0.490386734)

ko <- c(0.199400383, 0.216045179)

max_length <- max(length(ctrl), length(sleep), length(ko))

length(ctrl) <- max_length
length(sleep) <- max_length
length(ko) <- max_length

my_data <- data.frame( 
  group = rep(c("Control", "Sleep", "KO"), each = max_length),
  weight = c(ctrl,  sleep, ko)
)
my_data <- na.omit(my_data)

to_graph <- group_by(my_data, group) %>%
  summarise(
    count = n(),
    mean = mean(weight, na.rm = TRUE),
    sd = sd(weight, na.rm = TRUE)
  )
to_graph

t.test(ctrl, sleep, var.equal = TRUE)
t.test(ctrl, ko, var.equal = TRUE)
t.test(ko, sleep, var.equal = TRUE)

#Check Significance Marker matches expected
ggboxplot(my_data, x = "group", y = "weight", 
          color = "group", palette = c("#00AFBB", "#E7B800"),
          ylab = "Weight", xlab = "Groups")+
  geom_signif(comparisons = list(c("Control","Sleep")),map_signif_level = TRUE)

#Make Graph
ggplot(to_graph)+
  geom_bar( aes(x=group, y=mean, fill=group), stat="identity", alpha=0.7) +
  ggtitle("Hippocampus")+
  ylab("AKAP11 Fluorescence (Normalized to DAPI)")+
  xlab("")+
  theme(legend.position = "none")+
  geom_errorbar( aes(x=group, ymin=mean-sd, ymax=mean+sd), width=0.4, alpha=0.9, linewidth=0.5)+
  geom_point(aes(x = group, y = weight), my_data)+
  if(FALSE){geom_signif(
    data = to_graph,
    aes(xmin = 1, xmax = 2, annotations = "NS", y_position = 0.5),
    textsize = 3, vjust = -0.2,
    manual = TRUE) +
  geom_signif(
    data = to_graph,
    aes(xmin = 2, xmax = 3, annotations = "*", y_position = 0.80),
    textsize = 3, vjust = -0.2,
    manual = TRUE) +
  geom_signif(
    data = to_graph,
    aes(xmin = 1, xmax = 3, annotations = "*", y_position = 1.8),
    textsize = 3, vjust = -0.2,
    manual = TRUE) #+
  }  
  
