library(tidyverse)
install.packages("ggsignif")
library(ggsignif)
library(ggpubr)

#Set up Normalized Control and Sleep Data
ctrl <- c(0.843578645,
          0.77229271,
          
          1.148133236,
          1.235995409
          
)
sleep <- c(1.029244973,
           0.792766052,1.422377473,
           1.333309068)

my_data <- data.frame( 
  group = rep(c("Control", "Sleep"), each = 4),
  weight = c(ctrl,  sleep)
)

to_graph <- group_by(my_data, group) %>%
  summarise(
    count = n(),
    mean = mean(weight, na.rm = TRUE),
    sd = sd(weight, na.rm = TRUE)
  )
to_graph

t.test(ctrl, sleep, var.equal = TRUE)

#Check Significance Marker matches expected
ggboxplot(my_data, x = "group", y = "weight", 
          color = "group", palette = c("#00AFBB", "#E7B800"),
          ylab = "Weight", xlab = "Groups")+
  geom_signif(comparisons = list(c("Control","Sleep")),map_signif_level = TRUE)

#Make Graph
ggplot(to_graph)+
  geom_bar( aes(x=group, y=mean, fill=group), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=group, ymin=mean-sd, ymax=mean+sd), width=0.4, alpha=0.9, linewidth=0.5)+
  geom_point(aes(x = group, y = weight), my_data)+
  geom_signif(
    data = to_graph,
    aes(xmin = 1, xmax = 2, annotations = "NS", y_position = 1.5),
    textsize = 3, vjust = -0.2,
    manual = TRUE) +
  ggtitle("p62")+
  ylab("p62 (Normalized to Actin)")+
  xlab("")+
  theme(legend.position = "none")
