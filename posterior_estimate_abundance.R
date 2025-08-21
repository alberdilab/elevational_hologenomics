#elevation predictions without using beta parameter. try to see if the lowest and higher abundances with 90% confidence overlap


# Construct the gradient for prediction
pred_elevation <- constructGradient(m,
                                    focalVariable = "Elevation",
                                    non.focalVariables = 1,
                                    ngrid = gradientlength) %>%
  predict(m, Gradient = ., expected = TRUE)

# Convert predictions into a dataframe
pred_df <- do.call(rbind, pred_elevation) %>%
  as.data.frame() %>%
  mutate(elevation = rep(gradient, 1000)) %>%
  pivot_longer(-c(elevation), names_to = "genome", values_to = "value")

# Identify the lowest and highest elevation
elev_min <- min(pred_df$elevation)
elev_max <- max(pred_df$elevation)

# Compute 90% confidence intervals for the lowest and highest elevations
comparison_df <- pred_df %>%
  filter(elevation %in% c(elev_min, elev_max)) %>%
  group_by(elevation, genome) %>%
  summarise(
    Abundance_Low90 = quantile(value, 0.05),  # 5th percentile
    Abundance_High90 = quantile(value, 0.95), # 95th percentile
    Mean_Abundance = mean(value),             # Mean abundance
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = elevation, values_from = c(Abundance_Low90, Abundance_High90, Mean_Abundance))

# Rename columns for clarity
colnames(comparison_df) <- c("genome", 
                             "Low_Abundance_Low90", "High_Abundance_Low90",
                             "Low_Abundance_High90", "High_Abundance_High90",
                             "Low_Mean_Abundance", "High_Mean_Abundance")

# Check if confidence intervals overlap
comparison_df <- comparison_df %>%
  mutate(Overlap = ifelse(Low_Abundance_High90 >= High_Abundance_Low90, "Overlap", "No Overlap"))

# Convert data back to long format for plotting
plot_df <- comparison_df %>%
  pivot_longer(cols = c(Low_Mean_Abundance, High_Mean_Abundance), 
               names_to = "Elevation", 
               values_to = "Mean_Abundance") %>%
  mutate(Elevation = ifelse(Elevation == "Low_Mean_Abundance", "Lowest Elevation", "Highest Elevation"))

# Plot with overlapping MAGs highlighted
ggplot(plot_df, aes(x = genome, y = Mean_Abundance, color = Overlap, shape = Elevation)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Mean abundance points
  geom_errorbar(aes(ymin = ifelse(Elevation == "Lowest Elevation", Low_Abundance_Low90, High_Abundance_Low90), 
                    ymax = ifelse(Elevation == "Lowest Elevation", Low_Abundance_High90, High_Abundance_High90)), 
                width = 0.2, position = position_dodge(width = 0.5)) +  # Error bars
  scale_color_manual(values = c("No Overlap" = "black", "Overlap" = "orange")) +  # Highlight overlaps in orange
  labs(
    x = "MAG",
    y = "Predicted Abundance (90% CI)",
    color = "Overlap",
    shape = "Elevation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##elevation predictions without using beta parameter. try to see if the lowest and higher abundances with 90% confidence overlap (plot phylum and not MAG)

pred_df <- do.call(rbind, pred_elevation) %>%
  as.data.frame() %>%
  mutate(elevation = rep(gradient, 1000)) %>%
  pivot_longer(-c(elevation), names_to = "genome", values_to = "value")


pred_df_2<-pred_df%>%
  left_join(genome_metadata, by = join_by(genome == genome)) %>%
  select("elevation", "genome", "value", "phylum")

# Identify the lowest and highest elevation
elev_min <- min(pred_df_2$elevation)
elev_max <- max(pred_df_2$elevation)

# Compute 90% confidence intervals for the lowest and highest elevations
comparison_df_2 <- pred_df_2 %>%
  filter(elevation %in% c(elev_min, elev_max)) %>%
  group_by(elevation, phylum) %>%
  summarise(
    Abundance_Low90 = quantile(value, 0.05),  # 5th percentile
    Abundance_High90 = quantile(value, 0.95), # 95th percentile
    Mean_Abundance = mean(value),             # Mean abundance
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = elevation, values_from = c(Abundance_Low90, Abundance_High90, Mean_Abundance))

# Rename columns for clarity
colnames(comparison_df_2) <- c("phylum", 
                             "Low_Abundance_Low90", "High_Abundance_Low90",
                             "Low_Abundance_High90", "High_Abundance_High90",
                             "Low_Mean_Abundance", "High_Mean_Abundance")

# Check if confidence intervals overlap
comparison_df_2 <- comparison_df_2 %>%
  mutate(Overlap = ifelse(Low_Abundance_High90 >= High_Abundance_Low90, "Overlap", "No Overlap"))

# Convert data back to long format for plotting
plot_df_2 <- comparison_df_2 %>%
  pivot_longer(cols = c(Low_Mean_Abundance, High_Mean_Abundance), 
               names_to = "Elevation", 
               values_to = "Mean_Abundance") %>%
  mutate(Elevation = ifelse(Elevation == "Low_Mean_Abundance", "Lowest Elevation", "Highest Elevation"))

# Plot with overlapping MAGs highlighted
ggplot(plot_df_2, aes(x = phylum, y = Mean_Abundance, color = Overlap, shape = Elevation)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Mean abundance points
  geom_errorbar(aes(ymin = ifelse(Elevation == "Lowest Elevation", Low_Abundance_Low90, High_Abundance_Low90), 
                    ymax = ifelse(Elevation == "Lowest Elevation", Low_Abundance_High90, High_Abundance_High90)), 
                width = 0.2, position = position_dodge(width = 0.5)) +  # Error bars
  scale_color_manual(values = c("No Overlap" = "black", "Overlap" = "orange")) +  # Highlight overlaps in orange
  labs(
    x = "Phylum",
    y = "Predicted Abundance (90% CI)",
    color = "Overlap",
    shape = "Elevation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





##plotting only the mags that don't show overlapping abundances and plot the trends

pred_elevation <- constructGradient(m,
                                    focalVariable = "Elevation",
                                    non.focalVariables = 1,
                                    ngrid = gradientlength) %>%
  predict(m, Gradient = ., expected = TRUE) %>%
  do.call(rbind,.) %>%
  as.data.frame() %>%
  mutate(elevation=rep(gradient,1000)) %>%
  pivot_longer(-c(elevation), names_to = "genome", values_to = "value")

pred_elevation %>%
  left_join(genome_metadata, by = join_by(genome == genome)) %>%
  filter(genome %in% c("EHA01453_bin.112", "EHA01467_bin.143", "EHA01472_bin.16", "EHA01546_bin.13", "EHA01550_bin.77", "EHA01556_bin.110")) %>%
  select("elevation", "genome", "value", "phylum")%>%
  ggplot(., aes(x = elevation, y = value, fill = genome, colour=phylum)) +
  scale_color_manual(name="Phylum",
                     breaks=c("p__Desulfobacterota", "p__Pseudomonadota"),
                     labels=c("Desulfobacterota","Pseudomonadota"),
                     values=c("#008000","#ebe082")) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.2)+
  theme_classic() +
  xlab("Elevation (m)") +
  ylab("Predicted Abundance")

#talk about their phylogenetic relationhip and functional features




###overlaping abundances with posterior estimates (not possible without using beta) not working ####

post_estimates_Beta <- getPostEstimate(m, parName = "Beta")

beta_estimates <- data.frame(
  MAG = colnames(post_estimates_Beta$mean),
  Beta_Mean = post_estimates_Beta$mean[1, ],  
  Beta_Support = post_estimates_Beta$support[1, ],  
  Beta_SupportNeg = post_estimates_Beta$supportNeg[1, ]
)


E_low <- min(m$XData$Elevation, na.rm = TRUE)
E_high <- max(m$XData$Elevation, na.rm = TRUE)

abundance_estimates <- beta_estimates %>%
  mutate(
    Abundance_Low = Beta_Mean * E_low,
    Abundance_High = Beta_Mean * E_high
  )

abundance_estimates <- abundance_estimates %>%
  mutate(
    Low_Abundance_Low90 = Abundance_Low - 1.645 * sd(Abundance_Low),
    High_Abundance_High90 = Abundance_High + 1.645 * sd(Abundance_High)
  )

abundance_estimates <- abundance_estimates %>%
  mutate(Overlap = ifelse(High_Abundance_High90 >= Low_Abundance_Low90, "Overlap", "No Overlap"))

ggplot(abundance_estimates, aes(x = MAG, y = Abundance_High, color = Overlap)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Low_Abundance_Low90, ymax = High_Abundance_High90), width = 0.2) +
  scale_color_manual(values = c("Overlap" = "orange", "No Overlap" = "black")) +
  labs(
    title = "Abundance Overlap Between Lowest and Highest Elevation",
    x = "MAG",
    y = "Estimated Abundance",
    color = "Overlap"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#overlaping abundances of the raw data (not using this ones)

library(tidyverse)

# Extract raw abundance and corresponding elevation values
library(tidyverse)

# Extract raw abundance and corresponding elevation values
raw_abundance <- as.data.frame(m$Y) %>%
  mutate(Elevation = m$XData$Elevation) %>%
  pivot_longer(-Elevation, names_to = "genome", values_to = "Abundance")

# Define lowest and highest elevation
low_elevation <- min(m$XData$Elevation)  
high_elevation <- max(m$XData$Elevation)

# Filter for low and high elevation points
filtered_abundance <- raw_abundance %>%
  filter(Elevation %in% c(low_elevation, high_elevation)) %>%
  mutate(ElevationCategory = ifelse(Elevation == low_elevation, "Low Elevation", "High Elevation"))


# Compute 90% confidence intervals for the lowest and highest elevations
comparison_df <- filtered_abundance %>%
  filter(Elevation %in% c(low_elevation, high_elevation)) %>%
  group_by(Elevation, genome) %>%
  summarise(
    Abundance_Low90 = quantile(Abundance, 0.05),  # 5th percentile
    Abundance_High90 = quantile(Abundance, 0.95), # 95th percentile
    Mean_Abundance = mean(Abundance),             # Mean abundance
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Elevation, values_from = c(Abundance_Low90, Abundance_High90, Mean_Abundance))

# Rename columns for clarity
colnames(comparison_df) <- c("genome", 
                             "Low_Abundance_Low90", "High_Abundance_Low90",
                             "Low_Abundance_High90", "High_Abundance_High90",
                             "Low_Mean_Abundance", "High_Mean_Abundance")

# Check if confidence intervals overlap
comparison_df <- comparison_df %>%
  mutate(Overlap = ifelse(Low_Abundance_High90 >= High_Abundance_Low90, "Overlap", "No Overlap"))

# Convert data back to long format for plotting
plot_df <- comparison_df %>%
  pivot_longer(cols = c(Low_Mean_Abundance, High_Mean_Abundance), 
               names_to = "Elevation", 
               values_to = "Mean_Abundance") %>%
  mutate(Elevation = ifelse(Elevation == "Low_Mean_Abundance", "Lowest Elevation", "Highest Elevation"))

# Plot with overlapping MAGs highlighted
ggplot(plot_df, aes(x = genome, y = Mean_Abundance, color = Overlap, shape = Elevation)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Mean abundance points
  geom_errorbar(aes(ymin = ifelse(Elevation == "Lowest Elevation", Low_Abundance_Low90, High_Abundance_Low90), 
                    ymax = ifelse(Elevation == "Lowest Elevation", Low_Abundance_High90, High_Abundance_High90)), 
                width = 0.2, position = position_dodge(width = 0.5)) +  # Error bars
  scale_color_manual(values = c("No Overlap" = "black", "Overlap" = "orange")) +  # Highlight overlaps in orange
  labs(
    x = "MAG",
    y = "Observed Abundance (90% CI)",
    color = "Overlap",
    shape = "Elevation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#predicted function of the mags linked to elevation

# Aggregate bundle-level GIFTs into the compound level
#genome_counts_filt_filt <- tibble::rownames_to_column(genome_counts_filt_filt, var = "genome")
GIFTs_elements_filtered <- elements_table[rownames(elements_table) %in% genome_counts_filt$genome, ]
GIFTs_elements_filtered <- as.data.frame(GIFTs_elements_filtered) %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)


GIFTs_elements_filtered$genomes <- rownames(GIFTs_elements_filtered)

# Step 2: Add column names as a new row at the top
colnames_row <- as.character(colnames(df))  # Get column names as a vector
df <- rbind(colnames_row, df)

# Step 3: Update row names after adding the new row
rownames(df) <- c("Column_Names", rownames(df)[-1])



GIFTs_elements_filtered_2<-GIFTs_elements_filtered %>%
  filter(genomes %in% c("EHA01453_bin.112", "EHA01467_bin.143", "EHA01472_bin.16", "EHA01546_bin.13", "EHA01550_bin.77", "EHA01556_bin.110"))

# Remove second column
GIFTs_elements_filtered_2 <- GIFTs_elements_filtered_2[, -151]

# Get community-weighed average GIFTs per sample
GIFTs_elements_community_2 <- to.community(GIFTs_elements_filtered_2, genome_counts_filt %>% column_to_rownames(., "genome") %>% tss(), GIFT_db)

GIFTs_elements_community_2 %>%
  as.data.frame() %>%
  rownames_to_column(var="sample") %>%
  pivot_longer(!sample,names_to="trait",values_to="gift") %>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  mutate(functionid = substr(trait, 1, 3)) %>%
  mutate(trait = case_when(
    trait %in% GIFT_db$Code_element ~ GIFT_db$Element[match(trait, GIFT_db$Code_element)],
    TRUE ~ trait
  )) %>%
  mutate(functionid = case_when(
    functionid %in% GIFT_db$Code_function ~ GIFT_db$Function[match(functionid, GIFT_db$Code_function)],
    TRUE ~ functionid
  )) %>%
  mutate(trait=factor(trait,levels=unique(GIFT_db$Element))) %>%
  mutate(functionid=factor(functionid,levels=unique(GIFT_db$Function))) %>%
  ggplot(aes(x=sample,y=trait,fill=gift)) +
  geom_tile(colour="white", linewidth=0.2)+
  scale_fill_gradientn(colours=rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da")))+
  facet_grid(functionid ~ Elevation, scales="free",space="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Traits",x="Samples",fill="GIFT")






#when they are link to elevation
community_elements <- pred_elevation %>%
  filter(genome %in% c("EHA01453_bin.112", "EHA01467_bin.143", "EHA01472_bin.16", "EHA01546_bin.13", "EHA01550_bin.77", "EHA01556_bin.110")) %>%
  group_by(elevation, genome) %>%
  mutate(row_id = row_number()) %>%
  pivot_wider(names_from = genome, values_from = value) %>%
  ungroup() %>%
  group_split(row_id) %>%
  as.list() %>%
  lapply(., FUN = function(x){x %>%
      select(-row_id) %>%
      column_to_rownames(var = "elevation") %>%
      as.data.frame() %>%
      exp() %>%
      t() %>%
      tss() %>%
      to.community(elements_table,.,GIFT_db) %>%
      as.data.frame() %>%
      rownames_to_column(var="elevation")
  })

calculate_slope <- function(x) {
  lm_fit <- lm(unlist(x) ~ seq_along(unlist(x)))
  coef(lm_fit)[2]
}

element_predictions <- map_dfc(community_elements, function(mat) {
  mat %>%
    column_to_rownames(var = "elevation") %>%
    t() %>%
    as.data.frame() %>%
    rowwise() %>%
    mutate(slope = calculate_slope(c_across(everything()))) %>%
    select(slope) }) %>%
  t() %>%
  as.data.frame() %>%
  set_names(colnames(community_elements[[1]])[-1]) %>%
  rownames_to_column(var="iteration") %>%
  pivot_longer(!iteration, names_to="trait",values_to="value") %>%
  group_by(trait) %>%
  summarise(mean=mean(value),
            p1 = quantile(value, probs = 0.1),
            p9 = quantile(value, probs = 0.9),
            positive_support = sum(value > 0)/1000,
            negative_support = sum(value < 0)/1000) %>%
  arrange(-positive_support)

unique_funct_db<- GIFT_db[c(2,4,5,6)] %>% 
  distinct(Code_element, .keep_all = TRUE)

element_predictions %>%
  filter(mean >0) %>%
  arrange(-positive_support) %>%
  filter(positive_support>=0.9) %>%
  left_join(unique_funct_db, by = join_by(trait == Code_element))%>%
  arrange(Domain,Function)%>%
  paged_table()

element_predictions %>%
  filter(mean <0) %>%
  arrange(-negative_support) %>%
  filter(negative_support>=0.9) %>%
  left_join(unique_funct_db, by = join_by(trait == Code_element))%>%
  arrange(Domain,Function)%>%
  paged_table()

#Positively associated
positive <- element_predictions %>%
  filter(mean >0) %>%
  arrange(mean) %>%
  filter(positive_support>=0.9) %>%
  select(-negative_support) %>%
  rename(support=positive_support)

#Negatively associated
negative <- element_predictions %>%
  filter(mean <0) %>%
  arrange(mean) %>%
  filter(negative_support>=0.9) %>%
  select(-positive_support) %>%
  rename(support=negative_support)

all_elements <- bind_rows(positive,negative) %>%
  left_join(GIFT_db,by=join_by(trait==Code_element)) %>%
  mutate(trait=factor(trait,levels=c(rev(positive$trait),rev(negative$trait)))) %>%
  mutate(Code_function=factor(Code_function)) %>%
  mutate(element_legend=str_c(trait," - ",Element)) %>%
  mutate(function_legend=str_c(Code_function," - ",Function)) %>%
  select(trait,mean,p1,p9,element_legend,function_legend) %>% 
  unique()

gift_colors <- read_tsv("data/gift_colors.tsv") %>% 
  mutate(legend=str_c(Code_function," - ",Function))  %>% 
  filter(legend %in% all_elements$function_legend)

all_elements %>%
  ggplot(aes(x=mean, y=fct_reorder(element_legend, mean), xmin=p1, xmax=p9, color=function_legend)) +
  geom_point() +
  geom_errorbar() +
  xlim(c(-0.15,0.15)) +
  geom_vline(xintercept=0) +
  scale_color_manual(values = gift_colors$Color) +
  theme_minimal() +
  labs(x="Regression coefficient",y="Functional trait")