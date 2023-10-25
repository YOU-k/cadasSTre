### select region ###
seu_eye1_list[["dbit"]] <- readRDS("/data/to/your/path/seu_dbit_1.rds")
seu_eye1_list$dbit$anno <- seu_eye1_list$dbit@assays$RNA@counts["Pmel",]

### create df_input data frame specifying regions for analysis ###
data.frame(platform=rep(c("bmk","stomic","slideseq","dbit"),each=3),
           xmin=c(rep(500,3),rep(530,3),rep(450,3),rep(360,3)),
           xmax=c(rep(800,3),rep(830,3),rep(750,3),rep(660,3)),
           ymin=c(seq(350,600,100),seq(300,500,100),seq(350,600,100),seq(650,750,50)),
           ymax=c(seq(400,650,100),seq(350,550,100),seq(400,650,100),seq(700,800,50)),
           rep=rep(c(1:3),4)) -> df_input

### loop through the platforms and regions ###
for (i in  2:length(seu_eye1_list)){
  seu_tmp <- seu_eye1_list[[i]]
  xmin=df_input$xmin[(i-2)*3+1]
  xmax=df_input$xmax[(i-2)*3+1]
  ymin=df_input$ymin[(i-2)*3+1]
  ymax=df_input$ymax[(i-2)*3+1]
  # extract spatial data and generate plots for the region
  tmp_df <- data.frame(xcoord=seu_tmp@reductions$spatial@cell.embeddings[,1],
                       ycoord=seu_tmp@reductions$spatial@cell.embeddings[,2],
                       anno=seu_tmp@assays$RNA@counts["Pmel",])
  # create and save the plot
  png(paste0("/your/select_region/",names(seu_eye1_list)[i],"_1.png"),units = "in", res=300,width = 4,height = 4)
  p1 <- FeaturePlot(seu_tmp,reduction = "spatial",features = c("Pmel"),max.cutoff = "q99",min.cutoff = "q05",
                   raster = FALSE,order = TRUE,pt.size = ps) + 
    theme_minimal() +labs(title="") +
    theme( panel.background = element_blank(),  # remove the plot background
           panel.grid.major = element_blank(),  # remove major grid lines
           panel.grid.minor = element_blank(),
           axis.line = element_line(color = "black"),  # add x and y axis
           axis.title = element_text(size = 12) )+ # remove minor grid lines
    geom_rect(aes(xmin = xmin, xmax = xmax-50, ymin = ymin+20, ymax = ymax+50), fill = NA, col="black", alpha = 0.2) +
    scale_color_viridis_c(direction = -1,alpha = 0.8) +coord_fixed()+labs(x="",y="")
  print(p1)
  dev.off()
  png(paste0("/your/select_region/subset",names(seu_eye1_list)[i],"_1.png"),units = "in", res=300,width = 4,height = 4)
  p2 <- FeaturePlot(seu_tmp,reduction = "spatial",features = c("Pmel"),max.cutoff = "q99",min.cutoff = "q05",
                    raster = FALSE,order = TRUE,pt.size = 0.8) + 
    xlim(xmin-300,xmax+300) +
    ylim(ymin,ymax+150) +
    geom_vline(xintercept = xmin, linetype="dashed",  size=0.8) +
    geom_vline(xintercept = xmax, linetype="dashed",  size=0.8) +
    
    scale_color_viridis_c(direction = -1)+  theme_void()+ theme(legend.position = "none")+
    coord_fixed() 
  print(p2)
  dev.off()
  
}

# Define a function 'calc_integral' to calculate the integral of data
calc_integral <- function(x, y) {
  diff_x = diff(x)
  avg_y = (y[-1] + y[-length(y)]) / 2
  sum(diff_x * avg_y)
}

# Initialize an empty data frame 'agg_df_all'
agg_df_all <- data.frame()

# Iterate through each row of the 'df_input' data frame
for (i in 1:nrow(df_input)) {
  seu_tmp <- seu_eye1_list[[df_input$platform[i]]]
  
  # Extract 'xcoord' and 'ycoord' data from the spatial embeddings
  data.frame(
    xcoord = seu_tmp@reductions$spatial@cell.embeddings[, 1],
    ycoord = seu_tmp@reductions$spatial@cell.embeddings[, 2]
  ) -> df
  
  # Create 'x_group' and 'y_group' by rounding coordinates
  df <- df %>%
    mutate(
      x_group = floor(xcoord / 10) * 10,
      y_group = floor(ycoord / 10) * 10
    ) %>% 
    dplyr::filter(
      xcoord < df_input$xmax[i] & xcoord > df_input$xmin[i] &
        ycoord < df_input$ymax[i] & ycoord > df_input$ymin[i]
    )
  
  # Extract 'pmel' data for 'Pmel' feature
  df$pmel <- seu_tmp@assays$RNA@counts["Pmel", match(rownames(df), colnames(seu_tmp))]
  
  # Create an aggregated data frame 'agg_df' with the 'pmel' sums
  agg_df <- df %>%
    group_by(x_group) %>%
    summarise(pmel = sum(pmel, na.rm = TRUE))
  
  # Add 'platform' information to 'agg_df'
  agg_df$platform <- df_input$platform[i]
  
  # Calculate the integral for the data
  integral <- calc_integral(agg_df$x_group, agg_df$pmel)
  
  # Create a new column 'normalized' with normalized y-values
  agg_df <- agg_df %>%
    mutate(normalized = pmel / integral)
  
  # Plot the data with the normalized y-values
  ggplot(agg_df, aes(x = x_group, y = normalized)) +
    geom_line()
  
  # Adjust 'x_group' values to be relative to the minimum value
  agg_df$x_group <- agg_df$x_group - min(agg_df$x_group)
  
  # Add the 'rep' information
  agg_df$rep <- df_input$rep[i]
  
  # Append the aggregated data frame to 'agg_df_all'
  rbind(agg_df, agg_df_all) -> agg_df_all
}

# Convert 'rep' to a factor in the 'agg_df_all' data frame
agg_df_all$rep <- as.factor(agg_df_all$rep)

# Initialize 'average_df_all' and 'fwhm_df_all' data frames
average_df_all <- data.frame()
fwhm_df_all <- data.frame()

# Iterate through unique platforms in 'agg_df_all'
for (p in unique(agg_df_all$platform)) {
  # Subset 'agg_df_all' for the current platform
  agg_df_all[agg_df_all$platform == p,] -> agg_df_p
  
  # Create a line plot 'p1' for the platform
  p1 <- ggplot(agg_df_p, aes(x = x_group, y = pmel, color = rep)) +
    geom_line() 
  print(p1)
  
  # Calculate the peak values for each replicate
  peak_df <- agg_df_p %>%
    group_by(rep) %>%
    summarise(peak_x = x_group[which.max(pmel)])
  
  # Create 'fwhm_df' for the platform
  agg_df_p %>%
    group_by(rep) %>%
    summarise(fwhm = calculate_fwhm(x_group, pmel)) -> tmp_fwhm
  
  # Add the 'platform' information to 'tmp_fwhm' and append it to 'fwhm_df_all'
  tmp_fwhm$platform <- p
  rbind(tmp_fwhm, fwhm_df_all) -> fwhm_df_all
  
  # Join the peak values with the original data frame
  aligned_df <- agg_df_p %>%
    inner_join(peak_df, by = "rep") %>%
    mutate(x_group_aligned = x_group - peak_x)
  
  # Calculate the average values
  average_df <- aligned_df %>%
    group_by(x_group_aligned) %>%
    summarise(averaged = mean(pmel, na.rm = TRUE))
  
  # Create a line plot 'p2' for the average data
  p2 <- ggplot(average_df, aes(x = x_group_aligned, y = averaged)) +
    geom_smooth(span = 0.2, se = FALSE) +
    theme_minimal() 
  print(p2)
  
  # Calculate the integral for the average data
  integral <- calc_integral(average_df$x_group_aligned, average_df$averaged)
  
  # Normalize the y-values and append them to 'average_df_all'
  average_df <- average_df %>%
    mutate(normalized = averaged / integral)
  average_df$platform <- p
  rbind(average_df_all, average_df) -> average_df_all
}

# Save the resulting data frames
average_df_all -> average_df_all_1
fwhm_df_all -> fwhm_df_all_1
saveRDS(fwhm_df_all_1,"/your/result/path/fwhm_df_all_1.rds")
saveRDS(fwhm_df_all_2,"/your/result/path/fwhm_df_all_2.rds")
saveRDS(average_df_all_1,"/your/result/path/average_df_all_1.rds")
saveRDS(average_df_all_2,"/your/result/path/average_df_all_2.rds")

# Modify the 'platform' column in 'average_df_all' based on platform names
average_df_all <- average_df_all %>%
  mutate(platform = case_when(
    platform == "pixel" ~ "Pixel-seq",
    platform == "slide" ~ "Slide-seqV1.5",
    platform == "stomic" ~ "Stereo-seq",
    platform == "hdst" ~ "HDST",
    
    TRUE ~ platform
  ))

# Create a line plot using 'ggplot' for the 'slideseq' platform
ggplot(agg_df_all[agg_df_all$platform == "slideseq",], 
       aes(x = x_group, y = normalized, col = rep)) +
  #geom_line(aes(col = platform)) +
  geom_smooth(span = 0.2, se = FALSE, aes(col = rep), linewidth = 0.8)

# Create a line plot with aligned 'x_group' values in 'average_df_all'
ggplot(average_df_all, 
       aes(x = x_group_aligned, y = normalized)) +
  #geom_line(aes(col = platform)) +
  geom_smooth(span = 0.3, se = FALSE, aes(col = platform), linewidth = 0.8) +
  #scale_color_manual(values = platform_col) +
  labs(x = "Aligned Distance", y = "Normalized expression of Slc17a7", col = "Platform") +
  theme_minimal()

# Modify the 'platform' column in 'fwhm_df_all' based on platform names
fwhm_df_all %>% 
  mutate(platform = case_when(
    platform == "pixel" ~ "Pixel-seq",
    platform == "slide" ~ "Slide-seqV1.5",
    platform == "stomic" ~ "Stereo-seq",
    platform == "hdst" ~ "HDST",
    
    TRUE ~ platform
  )) %>% 
  ggplot()+ 
  # Create a boxplot and jitter plot for FWHM values
  geom_boxplot(aes(x = reorder(platform, fwhm), y = fwhm), outlier.colour = NA) +
  geom_jitter(aes(x = reorder(platform, fwhm), y = fwhm, col = platform)) +
  # Customize the appearance and labels
  theme_bw() +
  ylab("FWHM") + xlab("")