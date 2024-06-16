install.packages("tidyverse")
library(tidyverse)
install.packages("ggsignif")
library(ggsignif)
install.packages("ggpubr")
library(ggpubr)

#Set up Normalized Control and Sleep Data
No_RS_a <- c(0.067545305,
             0.081174439,
             0.032810271,
             0.041397154,
             0.034993271,
             0.034653465
             
             )

RS_a <- c(0.139303483,
          0.15625,
          0.045234249,
          0.037790698,
          0.057471264,
          0.066914498,
          0.048640916
          
)

No_RS_b <- c(0.06539075,
             0.074712644,
             0.033216783,
             0.039848197,
             0.063829787,
             0.045676998
             
             
)

RS_b <- c(0.05380334,
          0.077557756,
          0.079872204,
          0.059105431,
          0.112627986,
          0.141269841
          
          
)


max_length <- max(length(No_RS_a), length(RS_a),length(No_RS_b), length(RS_b))

length(No_RS_a) <- max_length
length(RS_a) <- max_length
length(No_RS_b) <- max_length
length(RS_b) <- max_length


my_data <- data.frame( 
  group = rep(c("Control A", "1 Hour RS A", "Control B", "1 Hour RS B"), each = max_length),
  weight = c(No_RS_a,  RS_a, No_RS_b, RS_b)
)
my_data <- na.omit(my_data)

to_graph <- group_by(my_data, group) %>%
  summarise(
    count = n(),
    mean = mean(weight, na.rm = TRUE),
    sd = sd(weight, na.rm = TRUE)
  )
to_graph

t.test(RS_a, No_RS_a, var.equal = TRUE)
t.test(RS_b, No_RS_b, var.equal = TRUE)
t.test(RS_a, RS_b, var.equal = TRUE)
t.test(No_RS_a, No_RS_b, var.equal = TRUE)

#Check Significance Marker matches expected
ggboxplot(my_data, x = "group", y = "weight", 
          color = "group", palette = c("#00AFBB", "#E7B800"),
          ylab = "Weight", xlab = "Groups")+
  geom_signif(comparisons = list(c("Control","Sleep")),map_signif_level = TRUE)

#Make Graph
ggplot(to_graph)+
  geom_bar( aes(x=group, y=mean, fill=group), stat="identity", alpha=0.7) +
  ggtitle("HIP")+
  ylab("CFOS/NeuN Count")+
  xlab("")+
  theme(legend.position = "none")+
  geom_errorbar( aes(x=group, ymin=mean-sd, ymax=mean+sd), width=0.4, alpha=0.9, linewidth=0.5)+
    geom_point(aes(x = group, y = weight), my_data)+
    geom_signif(
      data = to_graph,
      aes(xmin = 1, xmax = 2, annotations = "NS", y_position = 0.17),
      textsize = 3, vjust = -0.2,
      manual = TRUE)+
  geom_signif(
    data = to_graph,
    aes(xmin = 3, xmax = 4, annotations = "NS", y_position = 0.1),
    textsize = 3, vjust = -0.2,
    manual = TRUE)
