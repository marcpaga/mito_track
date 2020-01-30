# Contains the final plots 
require(tiff)
require(ggplot2)
require(RColorBrewer)

######## ROSE PLOTS ############################################################

# load the data for all the relevant movies and join it together
rose_coords <- read.csv(file = 'rose_coords.csv', header = TRUE)

temp_list <- list()

for (i in seq_len(nrow(rose_coords))) {
  print(i)
  video_name <- rose_coords$Rename[i]
  
  file_name <- paste0('out/tracking/', video_name, '.csv')
  plot_folder <- paste0('plots/', video_name, '/')
  tiff_file <- paste0('dat/', video_name, '.tif')
  dir.create(paste0("plots/", video_name))
  
  # read data
  raw_file <- read.csv(file = file_name, stringsAsFactors = FALSE)
  #img <- readTIFF(tiff_file, all = TRUE)
  
  # remove mito trackers that are not in all frames
  long_ids <- which(table(raw_file$id) == max(table(raw_file$id)))
  trace_df <- subset(raw_file, id %in% long_ids)
  
  #pixel sizes
  #x_coord <- c(1:dim(img[[1]])[2])
  #y_coord <- c(1:dim(img[[1]])[1])
  
  if ('dist_memb' %in% colnames(trace_df)) {
    trace_df$dist_memb <- NULL
  }
  
  trace_df$video_name <- video_name
  
  temp_list[[i]] <- trace_df
  
}

all_trace_df <- do.call('rbind', temp_list)

#### WT separate plots ####

for (vid in c('flow_wt_1', 'flow_wt_2', 'flow_wt_3', 'flow_wt_4')) {
  select_vid <- vid
  plotDf <- subset(all_trace_df, video_name == select_vid)
  plotDf <- plotDf[,2:13]
  
  target_area <- as.vector(rose_coords[rose_coords$Rename == select_vid, 2:5, drop = TRUE])
  
  selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]])
  selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])
  
  selected_trace$angle <- selected_trace$angle + 90
  selected_trace$angle[which(selected_trace$angle > 360)] <- selected_trace$angle[which(selected_trace$angle > 360)] - 360
  
  p <- ggplot(data = selected_trace, aes(x = angle, y = (..count..)/sum(..count..), fill = (..count..)/sum(..count..))) +
    geom_hline(yintercept = seq(0, 0.2, by = 0.05), colour = "grey30", size = 0.2) +
    geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "grey30", size = 0.2) +
    geom_histogram(breaks = seq(0, 360, 15), alpha = 0.8, color = 'black', size = 1) +
    coord_polar(start = 0) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, 'BuPu'),
                         limits = c(0, 0.2),
                         guide = guide_colorbar(title = 'Relative counts',
                                                frame.colour = 'black', 
                                                frame.linewidth = 2, 
                                                ticks.colour = 'black', 
                                                ticks.linewidth = 2,
                                                barheight = unit(4, 'cm'), 
                                                barwidth = unit(0.5, 'cm'), 
                                                title.theme = element_text(size = 12, face = 'bold'))) +
    scale_x_continuous(limits = c(0, 360), 
                       breaks = seq(0, 315, 45), 
                       labels = seq(0, 315, 45)) + 
    scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(colour = 'black',face = 'bold', size = 14), 
          legend.text = element_text(colour = 'black',face = 'bold', size = 12)) + 
    ggtitle(select_vid)
  
  
  ggsave(filename = paste0("plots/final/rose_", select_vid, '.pdf'), device = 'pdf')
  
}

#### WT joined plot ####

plotDf <- subset(all_trace_df, video_name %in% c('flow_wt_1', 'flow_wt_2', 'flow_wt_3', 'flow_wt_4'))
plotDf <- plotDf[,2:13]

target_area <- as.vector(rose_coords[rose_coords$Rename == select_vid, 2:5, drop = TRUE])

selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]])
selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])

selected_trace$angle <- selected_trace$angle + 90
selected_trace$angle[which(selected_trace$angle > 360)] <- selected_trace$angle[which(selected_trace$angle > 360)] - 360

p <- ggplot(data = selected_trace, aes(x = angle, y = (..count..)/sum(..count..), fill = (..count..)/sum(..count..))) +
  geom_hline(yintercept = seq(0, 0.2, by = 0.05), colour = "grey30", size = 0.2) +
  geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "grey30", size = 0.2) +
  geom_histogram(breaks = seq(0, 360, 15), alpha = 0.8, color = 'black', size = 1) +
  coord_polar(start = 0) +
  scale_fill_gradientn(colours = brewer.pal(n = 9, 'BuPu'),
                       limits = c(0, 0.2),
                       guide = guide_colorbar(title = 'Relative counts',
                                              frame.colour = 'black', 
                                              frame.linewidth = 2, 
                                              ticks.colour = 'black', 
                                              ticks.linewidth = 2,
                                              barheight = unit(4, 'cm'), 
                                              barwidth = unit(0.5, 'cm'), 
                                              title.theme = element_text(size = 12, face = 'bold'))) +
  scale_x_continuous(limits = c(0, 360), 
                     breaks = seq(0, 315, 45), 
                     labels = seq(0, 315, 45)) + 
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(colour = 'black',face = 'bold', size = 14), 
        legend.text = element_text(colour = 'black',face = 'bold', size = 12)) + 
  ggtitle('WT')


ggsave(filename = "plots/final/rose_wt.pdf", device = 'pdf')
  



#### NO bleach separate plots ####

for (vid in c('2_no_bleaching', '3_no_bleaching', '8_no_bleaching', '10_no_bleaching')) {
  select_vid <- vid
  plotDf <- subset(all_trace_df, video_name == select_vid)
  plotDf <- plotDf[,2:13]
  
  target_area <- as.vector(rose_coords[rose_coords$Rename == select_vid, 2:5, drop = TRUE])
  
  selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]])
  selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])
  
  selected_trace$angle <- selected_trace$angle + 90
  selected_trace$angle[which(selected_trace$angle > 360)] <- selected_trace$angle[which(selected_trace$angle > 360)] - 360
  
  p <- ggplot(data = selected_trace, aes(x = angle, y = (..count..)/sum(..count..), fill = (..count..)/sum(..count..))) +
    geom_hline(yintercept = seq(0, 0.2, by = 0.05), colour = "grey30", size = 0.2) +
    geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "grey30", size = 0.2) +
    geom_histogram(breaks = seq(0, 360, 15), alpha = 0.8, color = 'black', size = 1) +
    coord_polar(start = 0) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, 'BuPu'),
                         limits = c(0, 0.2),
                         guide = guide_colorbar(title = 'Relative counts',
                                                frame.colour = 'black', 
                                                frame.linewidth = 2, 
                                                ticks.colour = 'black', 
                                                ticks.linewidth = 2,
                                                barheight = unit(4, 'cm'), 
                                                barwidth = unit(0.5, 'cm'), 
                                                title.theme = element_text(size = 12, face = 'bold'))) +
    scale_x_continuous(limits = c(0, 360), 
                       breaks = seq(0, 315, 45), 
                       labels = seq(0, 315, 45)) + 
    scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(colour = 'black',face = 'bold', size = 14), 
          legend.text = element_text(colour = 'black',face = 'bold', size = 12)) + 
    ggtitle(select_vid)
  
  
  ggsave(filename = paste0("plots/final/rose_", select_vid, '.pdf'), device = 'pdf')
  
}


#### NO bleach joined plot ####

plotDf <- subset(all_trace_df, video_name %in% c('2_no_bleaching', '3_no_bleaching', '8_no_bleaching', '10_no_bleaching'))
plotDf <- plotDf[,2:13]

target_area <- as.vector(rose_coords[rose_coords$Rename == select_vid, 2:5, drop = TRUE])

selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]])
selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])

selected_trace$angle <- selected_trace$angle + 90
selected_trace$angle[which(selected_trace$angle > 360)] <- selected_trace$angle[which(selected_trace$angle > 360)] - 360

p <- ggplot(data = selected_trace, aes(x = angle, y = (..count..)/sum(..count..), fill = (..count..)/sum(..count..))) +
  geom_hline(yintercept = seq(0, 0.2, by = 0.05), colour = "grey30", size = 0.2) +
  geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "grey30", size = 0.2) +
  geom_histogram(breaks = seq(0, 360, 15), alpha = 0.8, color = 'black', size = 1) +
  coord_polar(start = 0) +
  scale_fill_gradientn(colours = brewer.pal(n = 9, 'BuPu'),
                       limits = c(0, 0.2),
                       guide = guide_colorbar(title = 'Relative counts',
                                              frame.colour = 'black', 
                                              frame.linewidth = 2, 
                                              ticks.colour = 'black', 
                                              ticks.linewidth = 2,
                                              barheight = unit(4, 'cm'), 
                                              barwidth = unit(0.5, 'cm'), 
                                              title.theme = element_text(size = 12, face = 'bold'))) +
  scale_x_continuous(limits = c(0, 360), 
                     breaks = seq(0, 315, 45), 
                     labels = seq(0, 315, 45)) + 
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(colour = 'black',face = 'bold', size = 14), 
        legend.text = element_text(colour = 'black',face = 'bold', size = 12)) + 
  ggtitle('No bleach')


ggsave(filename = "plots/final/rose_no_bleach.pdf", device = 'pdf')




#### bleach separate plots ####

for (vid in c('bleaching_1', 'bleaching_3', 'bleaching_6', 'bleaching_9')) {
  select_vid <- vid
  plotDf <- subset(all_trace_df, video_name == select_vid)
  plotDf <- plotDf[,2:13]
  
  target_area <- as.vector(rose_coords[rose_coords$Rename == select_vid, 2:5, drop = TRUE])
  
  selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]])
  selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])
  
  selected_trace$angle <- selected_trace$angle + 90
  selected_trace$angle[which(selected_trace$angle > 360)] <- selected_trace$angle[which(selected_trace$angle > 360)] - 360
  
  p <- ggplot(data = selected_trace, aes(x = angle, y = (..count..)/sum(..count..), fill = (..count..)/sum(..count..))) +
    geom_hline(yintercept = seq(0, 0.2, by = 0.05), colour = "grey30", size = 0.2) +
    geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "grey30", size = 0.2) +
    geom_histogram(breaks = seq(0, 360, 15), alpha = 0.8, color = 'black', size = 1) +
    coord_polar(start = 0) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, 'BuPu'),
                         limits = c(0, 0.2),
                         guide = guide_colorbar(title = 'Relative counts',
                                                frame.colour = 'black', 
                                                frame.linewidth = 2, 
                                                ticks.colour = 'black', 
                                                ticks.linewidth = 2,
                                                barheight = unit(4, 'cm'), 
                                                barwidth = unit(0.5, 'cm'), 
                                                title.theme = element_text(size = 12, face = 'bold'))) +
    scale_x_continuous(limits = c(0, 360), 
                       breaks = seq(0, 315, 45), 
                       labels = seq(0, 315, 45)) + 
    scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(colour = 'black',face = 'bold', size = 14), 
          legend.text = element_text(colour = 'black',face = 'bold', size = 12)) + 
    ggtitle(select_vid)
  
  
  ggsave(filename = paste0("plots/final/rose_", select_vid, '.pdf'), device = 'pdf')
  
}

#### Bleach joined plot ####

plotDf <- subset(all_trace_df, video_name %in% c('bleaching_1', 'bleaching_3', 'bleaching_6', 'bleaching_9'))
plotDf <- plotDf[,2:13]

target_area <- as.vector(rose_coords[rose_coords$Rename == select_vid, 2:5, drop = TRUE])

selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]])
selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])

selected_trace$angle <- selected_trace$angle + 90
selected_trace$angle[which(selected_trace$angle > 360)] <- selected_trace$angle[which(selected_trace$angle > 360)] - 360

p <- ggplot(data = selected_trace, aes(x = angle, y = (..count..)/sum(..count..), fill = (..count..)/sum(..count..))) +
  geom_hline(yintercept = seq(0, 0.2, by = 0.05), colour = "grey30", size = 0.2) +
  geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "grey30", size = 0.2) +
  geom_histogram(breaks = seq(0, 360, 15), alpha = 0.8, color = 'black', size = 1) +
  coord_polar(start = 0) +
  scale_fill_gradientn(colours = brewer.pal(n = 9, 'BuPu'),
                       limits = c(0, 0.2),
                       guide = guide_colorbar(title = 'Relative counts',
                                              frame.colour = 'black', 
                                              frame.linewidth = 2, 
                                              ticks.colour = 'black', 
                                              ticks.linewidth = 2,
                                              barheight = unit(4, 'cm'), 
                                              barwidth = unit(0.5, 'cm'), 
                                              title.theme = element_text(size = 12, face = 'bold'))) +
  scale_x_continuous(limits = c(0, 360), 
                     breaks = seq(0, 315, 45), 
                     labels = seq(0, 315, 45)) + 
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(colour = 'black',face = 'bold', size = 14), 
        legend.text = element_text(colour = 'black',face = 'bold', size = 12)) + 
  ggtitle('No bleach')


ggsave(filename = "plots/final/rose_bleach.pdf", device = 'pdf')


######## DISTANCE TO MEMBRANE PLOTS ############################################

memb_coords <- read.csv(file = 'memb_coords.csv', header = TRUE)

temp_list <- list()

for (i in seq_len(nrow(memb_coords))) {
  print(i)
  video_name <- memb_coords$Rename[i]
  
  file_name <- paste0('out/out/tracking/', video_name, '.csv')
  plot_folder <- paste0('plots/', video_name, '/')
  tiff_file <- paste0('dat/', video_name, '.tif')
  dir.create(paste0("plots/", video_name))
  
  # read data
  raw_file <- read.csv(file = file_name, stringsAsFactors = FALSE)
  #img <- readTIFF(tiff_file, all = TRUE)
  
  # remove mito trackers that are not in all frames
  long_ids <- which(table(raw_file$id) == max(table(raw_file$id)))
  trace_df <- subset(raw_file, id %in% long_ids)
  
  #pixel sizes
  #x_coord <- c(1:dim(img[[1]])[2])
  #y_coord <- c(1:dim(img[[1]])[1])

  
  trace_df$video_name <- video_name
  
  temp_list[[i]] <- trace_df
  
}

all_memb_df <- do.call('rbind', temp_list)


vids <- c(c('bleaching_1', 'bleaching_3', 'bleaching_6', 'bleaching_9'), 
          c('2_no_bleaching', '3_no_bleaching', '8_no_bleaching', '10_no_bleaching'))

colors <- c('Yes' = brewer.pal(11, 'PRGn')[2], 'No' = brewer.pal(11, 'PRGn')[10])
colors <- rep(colors, each = 4)

plot_list <- list()
for (i in seq_along(vids)) {
  select_vid <- vids[i]
  plotDf <- subset(all_memb_df, video_name == select_vid) 
  
  target_area <- as.vector(memb_coords[memb_coords$Rename == select_vid, 2:5, drop = TRUE])
  
  selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]] & t == 1)
  selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])
  
  plotDf <- subset(plotDf, id %in% selected_trace$id)
  
  p <- ggplot(data = plotDf) + 
    geom_path(aes(x = t, y = dist_memb, group = id), alpha = 0.05) +
    geom_smooth(aes(x = t, y = dist_memb), method = "lm", formula = y ~ splines::bs(x, 20), se = TRUE, color = colors[i]) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text = element_text(face = 'bold', size = 12, color = 'black'),
          axis.title = element_text(face = 'bold', size = 12, color = 'black')) + 
    coord_cartesian(xlim = c(min(all_memb_df$t), max(all_memb_df$t)),
                    ylim = c(0, 60)) + 
    labs(x = 'Time (s)', y = 'Distance to membrane (px)') + 
    ggtitle(select_vid)
  plot_list[[i]] <- p
}

cowplot::plot_grid(plotlist = plot_list, ncol = 2, nrow = 4)
ggsave(filename = "plots/final/dist_memb_separate.pdf", device = 'pdf')



plot_list <- list()
for (i in seq_along(vids)) {
  select_vid <- vids[i]
  plotDf <- subset(all_memb_df, video_name == select_vid) 
  
  target_area <- as.vector(memb_coords[memb_coords$Rename == select_vid, 2:5, drop = TRUE])
  
  selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]] & t == 1)
  selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])
  
  plotDf <- subset(plotDf, id %in% selected_trace$id)
  
  plotDf$id2 <- paste0(plotDf$id, '_', plotDf$video_name)
  plot_list[[i]] <- plotDf
}

plotDf <- do.call('rbind', plot_list)

plotDf$bleach <- NA
plotDf$bleach[which(grepl('no', plotDf$video_name))] <- 'No'
plotDf$bleach[which(!grepl('no', plotDf$video_name))] <- 'Yes'

p <- ggplot(data = plotDf) +
  geom_path(aes(x = t, y = dist_memb, group = id2, color = bleach), alpha = 0.1) +
  geom_smooth(aes(x = t, y = dist_memb, group = video_name, color = bleach),
              method = "lm", formula = y ~ splines::bs(x, 20), se = TRUE) +
  scale_color_manual(values = colors) + 
  theme_bw() +
  theme(axis.text = element_text(face = 'bold', size = 12, color = 'black'),
        axis.title = element_text(face = 'bold', size = 12, color = 'black')) + 
  coord_cartesian(xlim = c(min(all_memb_df$t), max(all_memb_df$t)),
                  ylim = c(0, 60)) + 
  labs(x = 'Time (s)', y = 'Distance to membrane (px)') 

ggsave(filename = "plots/final/dist_memb_together.pdf", device = 'pdf')





plot_list <- list()
for (i in seq_along(vids)) {
  select_vid <- vids[i]
  plotDf <- subset(all_memb_df, video_name == select_vid) 
  
  target_area <- as.vector(memb_coords[memb_coords$Rename == select_vid, 2:5, drop = TRUE])
  
  selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]] & t == 1)
  selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])
  
  plotDf <- subset(plotDf, id %in% selected_trace$id)
  
  plotDf$id2 <- paste0(plotDf$id, '_', plotDf$video_name)
  plot_list[[i]] <- plotDf
}

plotDf <- do.call('rbind', plot_list)

plotDf$bleach <- NA
plotDf$bleach[which(grepl('no', plotDf$video_name))] <- 'No'
plotDf$bleach[which(!grepl('no', plotDf$video_name))] <- 'Yes'

plotDf$t <- plotDf$t * 2
plotDf <- subset(plotDf, t < 301)

p <- ggplot(data = plotDf) +
  geom_path(aes(x = t, y = dist_memb, group = id2, color = bleach), alpha = 0.1) +
  geom_smooth(aes(x = t, y = dist_memb, group = bleach, color = bleach),
              method = "lm", formula = y ~ splines::bs(x, 20), se = TRUE) +
  scale_color_manual(values = colors, 
                     guide = guide_legend(title = 'Activation',
                                            frame.colour = 'black', 
                                            frame.linewidth = 2, 
                                            ticks.colour = 'black', 
                                            ticks.linewidth = 2,
                                            barheight = unit(4, 'cm'), 
                                            barwidth = unit(0.5, 'cm'), 
                                            title.theme = element_text(size = 12, face = 'bold'))) + 
  theme_bw() +
  theme(axis.text = element_text(face = 'bold', size = 12, color = 'black'),
        axis.title = element_text(face = 'bold', size = 12, color = 'black'),
        axis.line = element_line(size = 1, color = 'black'),
        axis.ticks = element_line(size = 1, color = 'black')) + 
  coord_cartesian(xlim = c(min(all_memb_df$t), 300),
                  ylim = c(0, 60), expand = c(0, 0)) + 
  labs(x = 'Time (s)', y = 'Distance to membrane (px)') 

ggsave(filename = "plots/final/dist_memb_together2.pdf", device = 'pdf')



# roseplot grid

as.character(rose_coords$Rename)
plot_list <- list()
for (vid in as.character(rose_coords$Rename)) {
  select_vid <- vid
  plotDf <- subset(all_trace_df, video_name == select_vid)
  plotDf <- plotDf[,2:13]
  
  target_area <- as.vector(rose_coords[rose_coords$Rename == select_vid, 2:5, drop = TRUE])
  
  selected_trace <- subset(plotDf, x > target_area[[1]] & x < target_area[[2]])
  selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])
  
  selected_trace$angle <- selected_trace$angle + 90
  selected_trace$angle[which(selected_trace$angle > 360)] <- selected_trace$angle[which(selected_trace$angle > 360)] - 360
  
  p <- ggplot(data = selected_trace, aes(x = angle, y = (..count..)/sum(..count..), fill = (..count..)/sum(..count..))) +
    geom_hline(yintercept = seq(0, 0.2, by = 0.05), colour = "grey30", size = 0.2) +
    geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "grey30", size = 0.2) +
    geom_histogram(breaks = seq(0, 360, 15), alpha = 0.8, color = 'black', size = 1) +
    coord_polar(start = 0) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, 'BuPu'),
                         limits = c(0, 0.2),
                         guide = guide_colorbar(title = 'Relative counts',
                                                frame.colour = 'black', 
                                                frame.linewidth = 2, 
                                                ticks.colour = 'black', 
                                                ticks.linewidth = 2,
                                                barheight = unit(4, 'cm'), 
                                                barwidth = unit(0.5, 'cm'), 
                                                title.theme = element_text(size = 12, face = 'bold'))) +
    scale_x_continuous(limits = c(0, 360), 
                       breaks = seq(0, 315, 45), 
                       labels = seq(0, 315, 45)) + 
    scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(colour = 'black',face = 'bold', size = 9), 
          legend.text = element_text(colour = 'black',face = 'bold', size = 12),
          legend.position = 'none') + 
    ggtitle(select_vid)
  
  
  plot_list[[vid]] <- p
  
}

plot_grid(plotlist = plot_list, ncol = 4, nrow = 3)
