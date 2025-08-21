#mag catalogue ring

# Generate  basal tree
circular_tree <- force.ultrametric(genome_tree, method="extend") %>% # extend to ultrametric for the sake of visualisation
  ggtree(., layout="fan", open.angle=10, size=0.5)

# Add phylum ring
circular_tree <- gheatmap(circular_tree, phylum_heatmap, offset=0, width=0.05, colnames=FALSE) +
  scale_fill_manual(values=phylum_colors) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), panel.margin = margin(0, 0, 0, 0)) +
  new_scale_fill()

# Add Aisa ring
circular_tree <- gheatmap(circular_tree, aisa_heatmap, offset=0.2, width=0.05, colnames=FALSE) +
  scale_fill_manual(values=c("#ffffff","#e5bd5b")) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), panel.margin = margin(0, 0, 0, 0)) +
  new_scale_fill()

# Add Aran ring
circular_tree <- gheatmap(circular_tree, aran_heatmap, offset=0.3, width=0.05, colnames=FALSE) +
  scale_fill_manual(values=c("#ffffff","#6b7398")) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), panel.margin = margin(0, 0, 0, 0)) +
  new_scale_fill()

# Add Sentein Verde ring
circular_tree <- gheatmap(circular_tree, sentein_heatmap, offset=0.4, width=0.05, colnames=FALSE) +
  scale_fill_manual(values=c("#ffffff","#e2815a")) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), panel.margin = margin(0, 0, 0, 0)) +
  new_scale_fill()

# Add Tourmalet ring
circular_tree <- gheatmap(circular_tree, tourmalet_heatmap, offset=0.5, width=0.05, colnames=FALSE) +
  scale_fill_manual(values=c("#ffffff","#876b96")) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), panel.margin = margin(0, 0, 0, 0)) +
  new_scale_fill()

# Add prevalence ring
circular_tree <-  circular_tree +
  new_scale_fill() +
  scale_fill_manual(values = "#cccccc") +
  geom_fruit(
    data=prevalence_data,
    geom=geom_bar,
    mapping = aes(x=prevalence, y=genome),
    offset = 0.4,
    orientation="y",
    stat="identity")

# Add completeness ring
circular_tree <- circular_tree +
  new_scale_fill() +
  scale_fill_gradient(low = "#d1f4ba", high = "#f4baba") +
  geom_fruit(
    data=genome_metadata,
    geom=geom_bar,
    mapping = aes(x=completeness, y=genome, fill=contamination),
    offset = 0.20,
    orientation="y",
    stat="identity")

# Add genome-size ring
circular_tree <-  circular_tree +
  new_scale_fill() +
  scale_fill_manual(values = "#cccccc") +
  geom_fruit(
    data=genome_metadata,
    geom=geom_bar,
    mapping = aes(x=length, y=genome),
    offset = 0.05,
    orientation="y",
    stat="identity")

# Add text
circular_tree <-  circular_tree +
  annotate('text', x=2.8, y=0, label='             Phylum', family='arial', size=3.5) +
  annotate('text', x=3.1, y=0, label='                Transect', family='arial', size=3.5) +
  annotate('text', x=3.8, y=0, label='                                Transect prevalence', family='arial', size=3.5) +
  annotate('text', x=5.0, y=0, label='                          Genome quality', family='arial',size=3.5)+
  annotate('text', x=5.4, y=0, label='                      Genome size', family='arial', size=3.5)

#Plot circular tree
circular_tree %>% open_tree(30) %>% rotate_tree(90)
