# Contains the final plots 
require(tiff)
require(ggplot2)
require(RColorBrewer)

video_name <- 'flow_wt_1'
# defines an area of interest 

# read the data
file_name <- paste0('data/processed/', video_name, '.csv')
plot_folder <- paste0('plots/', video_name, '/')
tiff_file <- paste0('data/rawdata/', video_name, '.tif')
dir.create(paste0("plots/", video_name))

raw_file <- read.csv(file = file_name, stringsAsFactors = FALSE)
img <- readTIFF(tiff_file, all = TRUE)

x_coord <- c(1:dim(img[[1]])[2])
y_coord <- c(1:dim(img[[1]])[1])



#### All mitochondria plot -----------------------------------------------------

# reverse y coordinates
trace_df <- raw_file
trace_df$y <- abs(trace_df$y - max(y_coord))

p <- ggplot() + 
  geom_path(data = trace_df, 
            aes(x = x, y = y, group = id, color = as.factor(id))) + 
  theme_bw() + 
  coord_cartesian(xlim = x_coord, ylim = y_coord) + 
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_blank())

ggsave(filename = paste0("plots/", video_name, '/global_movement.png'))


#### All mitochondria plot with targeted emphasis ------------------------------

target_area <- c(xmin = 300, xmax = 400, ymin = 100, ymax = 200)

trace_df <- raw_file
trace_df$y <- abs(trace_df$y - max(y_coord))
t1df <- subset(trace_df, t == 1)
t1df <- subset(t1df, x > target_area[1] & x < target_area[2])
t1df <- subset(t1df, y > target_area[3] & y < target_area[4])
selected_trace <- subset(trace_df, id %in% t1df$id)
non_selected_trace <- subset(trace_df, !(id %in% t1df$id))

p <- ggplot() + 
  geom_path(data = selected_trace, 
            aes(x = x, y = y, group = id, color = as.factor(id))) + 
  geom_path(data = non_selected_trace, 
            aes(x = x, y = y, group = id, color = as.factor(id)), alpha = 0.1) + 
  theme_bw() + 
  coord_cartesian(xlim = x_coord, ylim = y_coord) + 
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_blank())

ggsave(filename = paste0("plots/", video_name, '/global_movement_area.png'))

#### Distance to the membrane --------------------------------------------------

trace_df <- raw_file
selected_trace <- subset(trace_df, x > target_area[[1]] & x < target_area[[2]] & t == 1)
selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])

plotDf <- subset(trace_df, id %in% selected_trace$id)
# 1 frame = 2s
plotDf$t <- plotDf$t*2

p <- ggplot(data = plotDf) + 
  geom_path(aes(x = t, y = dist_memb, group = id), alpha = 0.05) +
  geom_smooth(aes(x = t, y = dist_memb), method = "lm", formula = y ~ splines::bs(x, 20), se = TRUE, color = 'black') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(face = 'bold', size = 12, color = 'black'),
        axis.title = element_text(face = 'bold', size = 12, color = 'black')) + 
  labs(x = 'Time (s)', y = 'Distance to membrane (px)')

ggsave(filename = paste0("plots/", video_name, '/distance_to_membrane_area.png'))

#### Rose plot -----------------------------------------------------------------


selected_trace <- subset(trace_df, x > target_area[[1]] & x < target_area[[2]])
selected_trace <- subset(selected_trace, y > target_area[[3]] & y < target_area[[4]])

## rotate angles so that 0 is at top
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
        legend.text = element_text(colour = 'black',face = 'bold', size = 12))

ggsave(filename = paste0("plots/", video_name, '/rose_plot_area.png'))