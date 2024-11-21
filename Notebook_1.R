# Fetch overall cell-level statistics (for all cell lines)
FDS_CT_stats <- readRDS("/data/AlzheimersDSData/overall_celltype_statistics.rds")
kable(head(FDS_CT_stats))

#Inspect the dimensions of this dataset using the base R dim() function:
dim(FDS_CT_stats)
## [1] 379142     3
#379142 rows, 3 columns
#In the data, we see that the columns are cell attributes:
#• CellType: The cell type category to which this cell belongs
#• Age: The developmental time point at which this cell was sequenced
#• Mt: The status of the cell as bearing the frontotemporal dementia disease mutation (V337M) or wildtype
#(V337V)


#Get counts of each cell type category:
t1 <- table(FDS_CT_stats$CellType)
kable(t1)

#Get counts of cells analyzed at each sequencing timepoint:
t2 <- table(FDS_CT_stats$Age)
kable(t2)

#Get counts at each combination of time point + CellType:
t3 <- table(FDS_CT_stats$CellType,FDS_CT_stats$Age)
kable(t3)

#As a visual alternative, we can plot the proportions of cells in each cell type, across the full dataset:
# Calculate overall cell type proportions for plotting
FDS_CT_props <- FDS_CT_stats %>%
  group_by(CellType,Age) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(Age) %>%
  mutate(freq = n / sum(n))
## `summarise()` has grouped output by 'CellType'. You can override using the
## `.groups` argument.
ggplot(FDS_CT_props, aes(x = Age, y=freq, fill = CellType)) +
  geom_col() +
  labs(title = "Cell Type Proportions at Each Timepoint",
       x = "Developmental Timepoint of Organoid Cells",
       y = "Proportion of Cells", fill="Cell Type" ) +
  theme_minimal()

#For each cell type, we can also inspect the distributions of cells with the 
#diseased phenotype vs. in the control 5 group
#(Recall: V337M = mutant, V337V = wildtype)

##Analyzing CellType-Specific Gene Expression Trends##
#Replace the code Ast (for “astrocytes”) below with the code for cell type that you would like to explore in
#depth. Make sure to keep it in double quotes (“ “)
# Important:
# The text you enter must exactly match
# the code for your cell type as given above- case-sensitive!
my_celltype <- "Ast"

#Replace the gene IDs below with the IDs of the marker genes for your cell type of interest
#Make sure to type the IDs exactly, keeping each one individually in double quotes (“ “).
# Important:
# The IDs you enter must exactly match the gene IDs for your cell type
# as given in the table- case-sensitive!
my_celltype_marker_genes <- c("GFAP","ALDH1L1","SOX9","AQP4")

#Read in the gene expression data and cell-level metadata for your cell type:
my_celltype_data <- readRDS(sprintf("/data/AlzheimersDSData/celltype_data_%s.rds", my_celltype))

#The function readRDS() reads the data from a specified file in your R environment into an R object. After
#read-in, the object can then be accessed for the remainder of the R session under the same variable name
# gather data into long format for plotting,
# with one column "Gene" containing all 3000 of the most-variable genes
# identified in the overall dataset
my_celltype_data.df <- my_celltype_data %>% tidyr::gather(key = "Gene",
                                                          value="ExpressionLevel",
                                                          -Age,-Mt,-CellType)
# filter the resulting data frame
# so that it contains only the rows corresponding to the
# cell type marker genes
my_ct_genes.df <- my_celltype_data.df %>%
  filter(Gene %in% my_celltype_marker_genes)
# Use `ggplot2` to plot expression distributions of marker genes for this cell type,
# by means of a "violin plot" (ggplot2::geom_violin())
plot <- ggplot(data = my_ct_genes.df,
               aes(x = Gene, y = ExpressionLevel, fill=Gene))
plot <- plot + geom_violin()
plot

#This visual is missing several key aesthetics that add value for data interpretation- most notably a title. We
#can use the ggplot2 library to add an informative title, and also change the plot theme to make it a bit
#more visually pleasing, simply by adding additional geoms (or “layers”) to our existing plot
plot <- plot + ggtitle(sprintf("Expression Distributions of Key Marker Genes for CellType %s", my_celltype)) +
                               theme_minimal()
plot

# select one marker gene for plotting
marker_gene1 <- my_celltype_marker_genes[3]
# filter the previous celltype data
# so that it contains only the rows corresponding to the
# cell type marker gene of interest
gene.df <- my_celltype_data.df %>%
  filter(Gene == marker_gene1)

#Plot expression distributions with respect to variant (attribute “Mt”)
var_to_plot <- "Mt"
# Use `ggplot` to plot the expression distribution of
# this cell type marker gene with respect to variant,
# by means of a "violin plot" (ggplot2::geom_violin())
plot2 <- ggplot(data = gene.df,
                aes(x = .data[[var_to_plot]], y = ExpressionLevel, fill=.data[[var_to_plot]]))
plot2 <- plot2 + geom_violin() +
  ggtitle(sprintf("%s Expression Distribution within CellType %s, Diseased vs. Normal Cell",my_celltype)) +
theme_minimal()
plot2

#Plot expression distributions with respect to time (attribute "Age")
var_to_plot <- "Age"
plot_age <- ggplot(data = gene.df, aes(x = .data[[var_to_plot]], y = ExpressionLevel, fill = .data[var_to_plot]))
plot_age <- plot_age + geom_violin + ggtitle(sprintf("%s Expression Distribution within Age %s, Diseased vs Normal Cells", my_celltype)) + theme_minimal()
plot_age

#Plot expression distributions with respect to CellType and time
var_to_plot <- "Age"
var_to_plot2 <- "Mt"
plot_age <- ggplot(data = gene.df, aes(x = .data[[var_to_plot]], y = ExpressionLevel, fill = .data[[var_to_plot2]]))
plot_age <- plot_age + geom_violin + ggtitle("%s Expression Distribution within CellType and Age %s, Diseased vs Normal Cells", my_celltype) + theme_minimal()
plot_age

#Alternative Visualization Methods
#render our expression distributions for the set of cell type marker genes as boxplots
# Use `ggplot` to plot expression distributions of marker genes for this cell type,
# this time by box plot (ggplot2::geom_boxplot())
plot3 <- ggplot(data = my_ct_genes.df,
                aes(x = Gene, y = ExpressionLevel, fill=Gene))
plot3 <- plot3 + geom_boxplot() +
  ggtitle(sprintf("Expression Distributions of Key Marker Genes for CellType %s",my_celltype)) + 
                  theme_minimal()
plot3

#We can even combine the aesthetics if desired!
plot + geom_boxplot(width=0.05, color="grey", alpha=0.6)

#An alternative option for studying the distribution of a numeric variable across several groups is a Ridgeplot,
#which we can similarly achieve by simply specifying a different ggplot geometry:

  # This graph is made using the ggridges library,
  # which is a ggplot2 extension
plot4 <- ggplot(data = my_ct_genes.df, aes(x = ExpressionLevel, y = Gene, fill=Gene))
plot4 <- plot4 + geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") + ggtitle(sprintf("Expression Distributions of Key Marker Genes for CellType Ast", my_celltype))
plot4

#, start with the data in my_ct_genes.df and come up with an alternative visual-
#ization geometry for rendering the expression distributions of your set of celltype marker genes
plot_dot <- ggplot(my_ct_genes.df, aes(x = ExpressionLevel, y = Gene))
plot_dot <- plot_dot + geom_point(aes(colour = Gene)) + theme_ridges() + ggtitle(sprintf("Expression Distributions of Key Marker Genes for CellType Ast", my_ct_genes.df)) + theme_ridges() + theme(legend.position = "none")
plot_dot
