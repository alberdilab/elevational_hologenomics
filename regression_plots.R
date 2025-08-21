#### Richness diversity

```{r alpha_div_richness_plot, message=FALSE, warning=FALSE, fig.height=4, fig.width=10, fig.fullwidth=TRUE}
columns <- c("richness","neutral","phylo","func","mapped","total")
alpha_div %>%
  select(sample,richness) %>%
  pivot_longer(-sample, names_to = "data", values_to = "value") %>%
  mutate(data = factor(data, levels = columns))	%>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  ggplot(aes(x = Elevation, y = value)) +
  geom_point() +
  geom_smooth(method = "loess") +
  stat_poly_eq()+
  facet_wrap(~ factor(Transect))+
  labs(x = "Elevation (m)")
```

#### Neutral diversity

```{r alpha_div_neutral_plot, message=FALSE, warning=FALSE, fig.height=4, fig.width=10, fig.fullwidth=TRUE}
alpha_div %>%
  select(sample,neutral) %>%
  pivot_longer(-sample, names_to = "data", values_to = "value") %>%
  mutate(data = factor(data, levels = columns))	%>%
  left_join(sample_metadata, by = join_by(sample == EHI_number))  %>%
  ggplot(aes(x = Elevation, y = value)) +
  geom_point() +
  geom_smooth(method = "loess") +
  stat_poly_eq()+
  facet_wrap(~ factor(Transect))+
  labs(x = "Elevation (m)")
```

#### Phylogenetic diversity

```{r alpha_div_phylo_plot, message=FALSE, warning=FALSE, fig.height=4, fig.width=10, fig.fullwidth=TRUE}
alpha_div %>%
  select(sample,phylogenetic) %>%
  pivot_longer(-sample, names_to = "data", values_to = "value") %>%
  mutate(data = factor(data, levels = columns))	%>%
  left_join(sample_metadata, by = join_by(sample == EHI_number)) %>%
  ggplot(aes(x = Elevation, y = value)) +
  geom_point() +
  geom_smooth(method = "loess") +
  stat_poly_eq()+
  facet_wrap(~ factor(Transect))+
  labs(x = "Elevation (m)")
```

#### Functional diversities

```{r alpha_div_func_plot, message=FALSE, warning=FALSE, fig.height=4, fig.width=10, fig.fullwidth=TRUE}
alpha_div %>%
  select(sample,functional) %>%
  pivot_longer(-sample, names_to = "data", values_to = "value") %>%
  mutate(data = factor(data, levels = columns))	%>%
  left_join(sample_metadata, by = join_by(sample == EHI_number))  %>%
  ggplot(aes(x = Elevation, y = value)) +
  geom_point() +
  geom_smooth(method = "loess") +
  stat_poly_eq()+
  facet_wrap(~ factor(Transect))+
  labs(x = "Elevation (m)")
```