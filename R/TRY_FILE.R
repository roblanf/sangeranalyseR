


dhc <- as.dendrogram(SangerAlignedConsensusSet@SCdendrogram[[2]])
# Build dendrogram object from hclust results
dend <- as.dendrogram(hc)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)


# Plot line segments and add labels
p <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    geom_text(data = dend_data$labels, aes(x, y, label = label),
              hjust = 1, angle = 0, size = 3) +
    coord_flip() +
    scale_y_reverse(expand = c(0.2, 0)) +
    theme_dendro()
ggplotly(p)

# Rectangular lines
p <- ggplot(segment(ddata)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = dend_data$labels, aes(x, y, label = label),
              hjust = 10, vjust = 10, angle = 0, size = 3, check_overlap = TRUE) +
    coord_flip() +
    scale_y_reverse(expand = c(0.2, 0)) +
    theme_dendro()

ggplotly(p)




pp <- ggdendrogram(dhc, rotate = TRUE)
ggplotly(pp)


hc <- hclust(dist(SangerAlignedConsensusSet@SCdendrogram[[1]]), "ave")
dend1 <- as.dendrogram(hc)

plot_dendro(dend1, height = "100%", width = "300%")


hc <- hclust(dist(A_chloroticConsensusReads@dendrogram[[1]]), "ave")
dend1 <- as.dendrogram(hc)
plot_dendro(dend1)


dend <- iris[1:30,-5] %>% scale %>% dist %>%
    hclust %>% as.dendrogram %>%
    set("branches_k_color", k=3) %>% set("branches_lwd", 1.2) %>%
    set("labels_colors") %>% set("labels_cex", c(.9,1.2)) %>%
    set("leaves_pch", 19) %>% set("leaves_col", c("blue", "red"))
# plot the dend in usual "base" plotting engine:
plot(dend)

