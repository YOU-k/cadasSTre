library(ggsci)

### create a color palette 'platform_col' with 7 colors using the 'nrc' palette from ggsci ###
pal_npg("nrc")(7) -> platform_col
### assign platform names to the 'platform_col' palette ###
names(platform_col) <- c("Dbit-seq","Visium","Stereo-seq","BMKMANU S1000","Slide-seq","Pixel-seq","HDST")
### create another color palette 'pf_col' with 7 colors using the 'nrc' palette from ggsci ###
pal_npg("nrc")(7) -> pf_col
### assign platform names to the 'pf_col' palette ###
names(pf_col) <- c("dbit","10x","stomic","bmk","slideseq","pixel","hdst")

### load your own dataset ###
seu_eye1_list <- readRDS("/path/to/your/data/seu_eye1_list.rds")
seu_eye2_list <- readRDS("/path/to/your/data/seu_eye2_list.rds")

### create data frames 'df1' and 'df2' with platform names, total counts, and replicate information ###
data.frame(platform=names(seu_eye1_list),
           total_counts=unlist(lapply(seu_eye1_list, function(seu){sum(seu$nCount_RNA)})),
           rep=1) -> df1
data.frame(platform=names(seu_eye2_list),
           total_counts=unlist(lapply(seu_eye2_list, function(seu){sum(seu$nCount_RNA)})),
           rep=2) -> df2  
### combine df1 and df2 into df_all ###
rbind(df1,df2) -> df_all

### create lists ###
seu_eye1_downsampled <- list()
seu_eye2_downsampled <- list()

### load and add data for different platform to the respective list ###
# 10x #
seu <- readRDS("/path/to/your/data/seu.rds")
seu_eye1_downsampled[["10x"]] <- seu
seu <- readRDS("/path/to/your/data/seu.rds")
seu_eye2_downsampled[["10x"]] <- seu
# bmk #
seu <- readRDS("/path/to/your/data/seu.rds")
seu_eye1_downsampled[["bmk"]] <- seu
seu <- readRDS("/path/to/your/data/seu.rds")
seu_eye2_downsampled[["bmk"]] <- seu

seu <- readRDS("/path/to/your/data/seu.rds")
seu_eye1_downsampled[["slideseq"]] <- seu

seu <- readRDS("/path/to/your/data/seu.rds")
seu_eye1_downsampled[["stomic"]] <- seu
seu <- readRDS("/path/to/your/data/seu.rds")
seu_eye2_downsampled[["stomic"]] <- seu

### save the downsampled data lists to RDS files ###
saveRDS(seu_eye1_downsampled,"/path/to/your/result/seu_eye1_downsampled.rds")
saveRDS(seu_eye2_downsampled,"/path/to/your/result/seu_eye2_downsampled.rds")

### create data frames 'df1' and 'df2' with platform names, total counts, and replicate information for downsampled data ###
data.frame(platform=names(seu_eye1_downsampled),
           total_counts=unlist(lapply(seu_eye1_downsampled, function(seu){sum(seu$nCount_RNA)})),
           rep=1) -> df1
data.frame(platform=names(seu_eye2_downsampled),
           total_counts=unlist(lapply(seu_eye2_downsampled, function(seu){sum(seu$nCount_RNA)})),
           rep=2) -> df2  

rbind(df1,df2) -> df_downsample
### create an 'id' column in 'df' by combining 'platform' and 'rep' columns ###
df$id <- paste0(df$platform, "_", df$rep)
### create an 'id' column in 'df_downsample' ###
df_downsample$id <- paste0(df_downsample$platform, "_", df_downsample$rep)
### rearrange the levels of 'df$id' based on total_counts in 'df_downsample' ###
df$id <- factor(df$id, levels = df_downsample$id[order(df_downsample$total_counts, decreasing = TRUE)])

### generate a PDF plot file for 'count_bar.pdf' ###
pdf("/mnt/disk2/yue_data/data/spatial_benchmark/eb/downsampled/plots/count_bar.pdf", width = 7, height = 2.5)

### plot a bar chart ###
df %>% mutate(pf = factor(df$platform, levels = c("dbit", "10x", "stomic", "bmk", "slideseq", "pixel", "hdst"),
                          labels = c("Dbit-seq", "Visium", "Stereo-seq", "BMKMANU S1000", "Slide-seq", "Pixel-seq", "HDST"))) %>% 
  ggplot(aes(x = id, y = total_counts, fill = pf)) + geom_col() + 
  facet_wrap(. ~ counts, nrow = 1) +  
  scale_fill_manual(values = platform_col) + theme_classic() +
  theme(
    axis.text.x = element_blank(),  # Set x-axis labels to be transparent to make them invisible
    axis.title.x = element_blank(),  # Set x-axis title to be transparent to make it invisible
    strip.background = element_blank(),  # Remove the background of facet labels
    # axis.line.x = element_blank(),  # Remove x-axis line
    axis.ticks.x = element_blank(), 
    strip.text = element_text(size = 12)
  ) +
  labs(x = "", y = "Total counts", fill = "Platforms")
### save the plot to the PDF file ###
dev.off()
### save the 'df' data frame to an RDS file ###
saveRDS(df, "/mnt/disk2/yue_data/data/spatial_benchmark/eb/downsampled/count_bar.rds")