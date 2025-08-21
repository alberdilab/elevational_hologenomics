# Beta diversity per elevation chunks

##Aisa transect

Aisa_transect_A<-Aisa_transect %>%
  filter(Elevation=="1250" | Elevation=="1450")

Aisa_counts_A <- temp_genome_counts[, which(colnames(temp_genome_counts) %in% rownames(Aisa_transect_A))]
identical(sort(colnames(Aisa_counts_A)), sort(as.character(rownames(Aisa_transect_A))))
Aisa_counts_func_A <- temp_genome_counts_func[, which(colnames(temp_genome_counts_func) %in% rownames(Aisa_transect_A))]

beta_div_richness_aisa_a<-hillpair(data=Aisa_counts_A, q=0)
beta_div_neutral_aisa_a<-hillpair(data=Aisa_counts_A, q=1)
beta_div_phylo_aisa_a<-hillpair(data=Aisa_counts_A, q=1, tree=genome_tree)
beta_div_func_aisa_a<-hillpair(data=Aisa_counts_func_A, q=1, dist=dist)

pairwise.adonis(beta_div_richness_aisa_a$S, Aisa_transect_A$Elevation, perm = 999)
pairwise.adonis(beta_div_neutral_aisa_a$S, Aisa_transect_A$Elevation, perm = 999)
pairwise.adonis(beta_div_phylo_aisa_a$S, Aisa_transect_A$Elevation, perm = 999)
pairwise.adonis(beta_div_func_aisa_a$S, Aisa_transect_A$Elevation, perm = 999)

##using population level beta div

### calculating all the chunks in the same go

Aisa_transect<-tibble::rownames_to_column(Aisa_transect, var = "sample")
Aisa_counts<-tibble::rownames_to_column(Aisa_counts, var = "genome")


richness_beta_aisa <- Aisa_counts %>%
  pivot_longer(!genome,names_to="sample",values_to = "counts") %>%
  left_join(Aisa_transect, by = join_by(sample == sample)) %>%
  mutate(Elevation=factor(Elevation)) %>% 
  group_by(Elevation) %>% 
  group_split() %>% 
  map_dbl(., ~ .x %>%
            select(genome, sample, counts) %>%
            pivot_wider(names_from = sample, values_from = counts) %>%
            column_to_rownames(var = "genome") %>%
            tss() %>%
            as.data.frame() %>%
            hilldiss(data=., metric="S", q = 0)
  )

richness_beta_aisa <- Aisa_counts %>%
  pivot_longer(!genome, names_to = "sample", values_to = "counts") %>%
  left_join(Aisa_transect, by = join_by(sample == sample)) %>%
  mutate(Elevation = factor(Elevation)) %>%
  group_by(Elevation) %>%
  group_split() %>%
  map(., ~ .x %>%
        select(genome, sample, counts) %>%
        pivot_wider(names_from = sample, values_from = counts) %>%
        column_to_rownames(var = "genome") %>%
        tss() %>%
        as.data.frame() %>%
        hilldiss(data = ., metric = "S", q = 0)) %>%
  do.call(cbind, .)


names(richness_beta_aisa)  <- unique(Aisa_transect$Elevation) %>% sort()

neutral_beta_aisa <- Aisa_counts %>%
  pivot_longer(!genome,names_to="sample",values_to = "counts") %>%
  left_join(Aisa_transect, by = join_by(sample == sample)) %>%
  mutate(Elevation=factor(Elevation)) %>% 
  group_by(Elevation) %>% 
  group_split() %>% 
  map_dbl(., ~ .x %>%
            select(genome, sample, counts) %>%
            pivot_wider(names_from = sample, values_from = counts) %>%
            column_to_rownames(var = "genome") %>%
            tss() %>%
            as.data.frame() %>%
            hilldiss(data=., metric="S", q = 1)
  )

names(neutral_beta_aisa)  <- unique(Aisa_transect$Elevation) %>% sort()

phylo_beta_aisa <- Aisa_counts %>%
  pivot_longer(!genome,names_to="sample",values_to = "counts") %>%
  left_join(Aisa_transect, by = join_by(sample == sample)) %>%
  mutate(Elevation=factor(Elevation)) %>% 
  group_by(Elevation) %>% 
  group_split() %>% 
  map_dbl(., ~ .x %>%
            select(genome, sample, counts) %>%
            pivot_wider(names_from = sample, values_from = counts) %>%
            column_to_rownames(var = "genome") %>%
            tss() %>%
            as.data.frame() %>%
            hilldiss(data=., metric="S", q = 1, tree=genome_tree)
  )

names(phylo_beta_aisa)  <- unique(Aisa_transect$Elevation) %>% sort()

Aisa_counts_func<-tibble::rownames_to_column(Aisa_counts_func, var = "genome")

func_beta_aisa <- Aisa_counts_func %>%
  pivot_longer(!genome,names_to="sample",values_to = "counts") %>%
  left_join(Aisa_transect, by = join_by(sample == sample)) %>%
  mutate(Elevation=factor(Elevation)) %>% 
  group_by(Elevation) %>% 
  group_split() %>% 
  map_dbl(., ~ .x %>%
            select(genome, sample, counts) %>%
            pivot_wider(names_from = sample, values_from = counts) %>%
            column_to_rownames(var = "genome") %>%
            tss() %>%
            as.data.frame() %>%
            hilldiss(data=., metric="S", q = 1, dist=dist)
  )

names(func_beta_aisa)  <- unique(Aisa_transect$Elevation) %>% sort()


# Merge all metrics
beta_div_aisa <- bind_rows(richness_beta_aisa,neutral_beta_aisa,phylo_beta_aisa,func_beta_aisa) %>%
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="Elevation") %>% 
  as_tibble() %>% 
  rename("richness"=V1,"neutral"=V2,"phylogenetic"=V3,"functional"=V4)


pairwise.adonis(beta_div_aisa$richness, beta_div_aisa$Elevation, perm = 999)

pairwise_comparisons(
  data = beta_div_aisa,
  x = richness,
  y = Elevation,
  type = "parametric",
  var.equal = TRUE,
  paired = FALSE,
  #p.adjust.method = "bonferroni"
)

###dividing the chunks


Aisa_transect_A<-tibble::rownames_to_column(Aisa_transect_A, var = "sample")

richness_beta_aisa_a <- genome_counts_filt %>%
  pivot_longer(!genome,names_to="sample",values_to = "counts") %>%
  left_join(Aisa_transect_A, by = join_by(sample == sample)) %>%
  mutate(Elevation=factor(Elevation)) %>% 
  group_by(Elevation) %>% 
  group_split() %>% 
  map_dbl(., ~ .x %>%
            select(genome, sample, counts) %>%
            pivot_wider(names_from = sample, values_from = counts) %>%
            column_to_rownames(var = "genome") %>%
            tss() %>%
            as.data.frame() %>%
            hilldiss(data=., metric="S", q = 0)
  )

pairwise.adonis(richness_beta_aisa_a, Aisa_transect_A$Elevation, perm = 999)



Aisa_transect_B<-Aisa_transect %>%
  filter(Elevation=="1450" | Elevation=="1650")

Aisa_counts_B <- temp_genome_counts[, which(colnames(temp_genome_counts) %in% rownames(Aisa_transect_B))]
identical(sort(colnames(Aisa_counts_B)), sort(as.character(rownames(Aisa_transect_B))))
Aisa_counts_func_B <- temp_genome_counts_func[, which(colnames(temp_genome_counts_func) %in% rownames(Aisa_transect_B))]

beta_div_richness_aisa<-hillpair(data=Aisa_counts, q=0)
beta_div_neutral_aisa<-hillpair(data=Aisa_counts, q=1)
beta_div_phylo_aisa<-hillpair(data=Aisa_counts, q=1, tree=genome_tree)
beta_div_func_aisa<-hillpair(data=Aisa_counts_func, q=1, dist=dist)

Aisa_transect_C<-Aisa_transect %>%
  filter(Elevation=="1650" | Elevation=="1850")

Aisa_counts_C <- temp_genome_counts[, which(colnames(temp_genome_counts) %in% rownames(Aisa_transect_C))]
identical(sort(colnames(Aisa_counts_C)), sort(as.character(rownames(Aisa_transect_C))))
Aisa_counts_func_C <- temp_genome_counts_func[, which(colnames(temp_genome_counts_func) %in% rownames(Aisa_transect_C))]

beta_div_richness_aisa<-hillpair(data=Aisa_counts, q=0)
beta_div_neutral_aisa<-hillpair(data=Aisa_counts, q=1)
beta_div_phylo_aisa<-hillpair(data=Aisa_counts, q=1, tree=genome_tree)
beta_div_func_aisa<-hillpair(data=Aisa_counts_func, q=1, dist=dist)



#alpha diversity correlations per transect

neutral_cor<-alpha_div %>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  filter(Transect=="Tourmalet") %>%
  lm(neutral ~ Elevation + I(Elevation^2)) %>%
  broom.mixed::tidy()


neutral_cor<-alpha_div %>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  filter(Transect=="Tourmalet") %>%
  lm(neutral ~ I(Elevation^2), data = ., REML = FALSE) %>%
  broom.mixed::tidy() %>%
  tt()

alpha_div %>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  filter(Transect=="Tourmalet") %>%
  ggplot(aes(Elevation, neutral), aes(x = Elevation, y = neutral)) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = "blue") + 
  stat_poly_eq()+
  labs(title = "Quadratic Regression: Alpha Diversity vs. Elevation", x = "Elevation", y = "Alpha Diversity")
  


model<-lm(neutral_cor$neutral ~ neutral_cor$Elevation + I(neutral_cor$Elevation^2))
  
summary(model)

model_gam <- gam(alpha_diversity ~ s(elevation))

model_poly <- lmerTest::lmer(neutral_cor$neutral ~ poly(neutral_cor$Elevation, 3)+(1 | neutral_cor$Sampling_point)) # Polynomial degree 2 # Summarize the model summary(model_poly)


alpha_div %>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  filter(Transect=="Tourmalet") %>%
  lm(alpha_diversity ~ poly(Elevation, 3))%>%
  summary()


richness_model<-alpha_div %>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  gam(richness ~ log(sequencing_depth) + s(Elevation, bs = "tp", k = 5) * Transect, data = .) %>%
  broom.mixed::tidy()

model_poly <- lm(richness ~ log(sequencing_depth) + poly(Elevation, 3) + Transect, data = your_data)


#plot elevation per population an elevations

beta_div %>%
  pivot_longer(-Sampling_point, names_to = "metric", values_to = "value") %>%
  left_join(sample_metadata %>% select(Transect,Sampling_point, Elevation) %>% unique, by = "Sampling_point") %>%
  filter(metric=="neutral") %>%
  ggplot(aes(y = value, x = Elevation ,group=Elevation, colour=Transect)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Transect), alpha = 0.5, size = 2) + 
  geom_line(aes(group = interaction(Transect, metric)), position = position_dodge(width = 0.5))+
  #scale_color_viridis_c(name = "Elevation", option = "viridis") +
  scale_color_manual(name="Transect",
                     breaks=c("Aisa","Aran","Sentein","Tourmalet"),
                     labels=c("Aisa","Aran","Sentein","Tourmalet"),
                     values=c("#e5bd5b50", "#6b739850","#e2815a50", "#876b9650")) +
  facet_wrap(. ~ Transect, scales = "free", ncol=4) +
  coord_cartesian(xlim = c(1, NA)) +
  #stat_compare_means(size=2) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + 
  ylab("Beta diversity")
