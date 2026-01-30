library(tidyverse)
library(readr)
library(broom)
library(cowplot)
library(caret)
library(tibble)
library(GGally)
library(stringr)
library(RColorBrewer)
library(forcats)
library(ggpubr)
library(lubridate)
library(reshape2)
library(xts) #interactive plots
library(RVAideMemoire)
library(gridExtra)
library(grid)
library(PMCMRplus)
library(car)
library(ggridges)
## add plotting libraries later

setwd("/work/hs325/env876")
df <- read.csv("/work/hs325/env876/raw.csv")
colnames(df) <- c("time", "temR1", "temR2", "temR3", "temR4", "temR5", "temR6")
head(df)

df_labeled <- df %>%
  mutate(
    time_clean = as.POSIXct((time - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
  ) %>%
  dplyr::rename(Deep = temR1, Shallow = temR6)

df_long <- df_labeled %>%
  select(time_clean, Deep, Shallow) %>%
  pivot_longer(cols = c(Deep, Shallow), names_to = "Depth", values_to = "Temp")

# 10 min subsample
df_long_10min <- df_long %>%
  dplyr::mutate(time_10min = lubridate::floor_date(time_clean, unit = "10 minutes")) %>%
  dplyr::group_by(time_10min, Depth) %>%
  dplyr::summarize(Temp = mean(Temp, na.rm = TRUE), .groups = "drop")

# hourly subsample
df_long_hour <- df_long %>%
  dplyr::mutate(time_hour = floor_date(time_clean, unit = "hour")) %>%
  dplyr::group_by(time_hour, Depth) %>%
  dplyr::summarize(Temp = mean(Temp, na.rm = TRUE), .groups = "drop")

# daily subsample
df_long_day <- df_long %>%
  dplyr::mutate(time_day = floor_date(time_clean, unit = "day")) %>%
  dplyr::group_by(time_day, Depth) %>%
  dplyr::summarize(Temp = mean(Temp, na.rm = TRUE), .groups = "drop")


df_long_10min$Depth <- factor(df_long_10min$Depth, levels = c("Shallow", "Deep"))
df_long_hour$Depth <- factor(df_long_hour$Depth, levels = c("Shallow", "Deep"))
df_long_day$Depth <- factor(df_long_day$Depth, levels = c("Shallow", "Deep"))


################################################################################

## Summary statistics
shallow <- subset(df_long_10min, Depth == "Shallow")
deep <- subset(df_long_10min, Depth == "Deep")

mean(shallow$Temp) # 27.94962 
range(shallow$Temp) # 19.87741 - 33.99100
sd(shallow$Temp) # 2.558511

mean(deep$Temp) # 27.82802
range(deep$Temp) # 20.80708 - 32.50217
sd(deep$Temp) # 2.450915

#KS test
ks.test(shallow$Temp,deep$Temp) 
#D = 0.03176, p-value < 2.2e-16

#Permanova on temperature distributions
perm.anova(Temp ~ Depth, data = df_long_10min, nperm = 999)
# Sum Sq     Df Mean Sq F value Pr(>F)    
# Depth         877      1  877.06  139.74  0.001 ***
#   Residuals 1489246 237274    6.28                

#Check equality of variances
leveneTest(Temp ~ Depth, data = df_long_10min, na.action = na.pass)
# Levene's Test for Homogeneity of Variance (center = median: na.pass)
#           Df F value    Pr(>F)    
# group      1  267.87 < 2.2e-16 ***
#       237274                      

#Percentage of time above 30C
# https://floridakeys.noaa.gov/coral-bleaching-faqs.html
df_long_10min %>%
  dplyr::mutate(exceed = Temp > 30) %>%
  dplyr::group_by(Depth) %>%
  dplyr::summarize(pct_time = mean(exceed) * 100)
# Deep 27.7
# Shallow 29.9

# Percentage of variance in May-Aug
df_stats <- df_long_10min %>%
  mutate(month = month(time_10min),
         is_summer = month %in% c(5, 6, 7, 8))
total_var <- var(df_stats$Temp, na.rm = TRUE)
summer_var <- var(df_stats$Temp[df_stats$is_summer], na.rm = TRUE)

# 3. Calculate percentage
percent_variability <- (summer_var / total_var) * 100
percent_variability
## 30.746% in May-July, 31.47442% in May-August

###################### Density plot

p1 <- ggdensity(df_long_10min, x = "Temp",
                add = "mean", rug = TRUE,
                color = "Depth", fill = "Depth",
                xlab = "Temperature (°C)", ylab = "Density",
                palette = c("#E7B800","#00AFBB")) +
  annotate("text", 
           x = min(df_long_10min$Temp, na.rm = TRUE), 
           y = Inf, 
           label = "--- Deep Mean", 
           color = "#00AFBB", 
           vjust = 4,      
           hjust = -0.1,   
           size = 3) +
  annotate("text", 
           x = min(df_long_10min$Temp, na.rm = TRUE), 
           y = Inf, 
           label = "--- Shallow Mean", 
           color = "#E7B800",
           vjust = 2,    
           hjust = -0.08,   
           size = 3) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 13), 
        axis.text = element_text(size = 11))
p1

###################### Line plot

p2 <- ggplot(df_long_10min, aes(x = time_10min, y = Temp, color = Depth)) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("Deep" = "#00AFBB", "Shallow" = "#E7B800")) +
  theme_cowplot() +
  annotate("text", 
           x = min(df_long_10min$time_10min, na.rm = TRUE), 
           y = -Inf,                                       
           label = "Two-Sided KS p < 0.001", 
           hjust = 0.05,                                  
           vjust = -1,                                    # Nudge up
           size = 3.5, 
           fontface = "italic") +
  labs(x = "Date", 
       y = "Temperature (°C)",
       color = "Sensor Depth") +
  theme(legend.position = "none", axis.title = element_text(size = 13), axis.text = element_text(size=11))
p2
#Two-Sided KS D = 0.03176, p < 0.001


###################### Corr plot

df_10min_wide <- df_long_10min %>%
  pivot_wider(names_from = Depth, values_from = Temp)
df_hour_wide <- df_long_hour %>%
  pivot_wider(names_from = Depth, values_from = Temp)
df_day_wide <- df_long_day %>%
  pivot_wider(names_from = Depth, values_from = Temp)

p3 <- ggplot(df_10min_wide, aes(x = Shallow, y = Deep)) +
  geom_point(alpha = 0.2, color = "steelblue") +
  geom_smooth(method = "lm", color = "firebrick", size = 2, se = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "spearman", 
           label.x.npc = "left", label.y.npc = "top") +
  labs(x = "Shallow Temperature (°C)",
       y = "Deep Temperature (°C)") +
  theme_half_open() + theme(axis.title = element_text(size = 13), axis.text = element_text(size=11))
p3

###### KS test on range of daily variation
range <- df_long_10min %>%
  dplyr::mutate(date = as.Date(time_10min)) %>%
  dplyr::group_by(date, Depth) %>%
  dplyr::summarize(
    daily_range = max(Temp) - min(Temp),
    .groups = "drop"
  )
shallow_range <- subset(range, Depth == "Shallow")
deep_range <- subset(range, Depth == "Deep")
ks.test(shallow_range$daily_range,deep_range$daily_range) 
# D = 0.50905, p-value < 2.2e-16


################ Line plot with daily ranges
p4 <- ggplot(range, aes(x = date, y = daily_range, color = Depth)) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("Deep" = "#00AFBB", "Shallow" = "#E7B800")) +
  theme_cowplot() +
  labs(x = "Date", 
       y = "Daily Temperature Range (°C)",
       color = "Sensor Depth") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size=11))
p4


################ Hysteresis plot
df_plot <- df_day_wide %>%
  arrange(time_day) %>%
  mutate(time_num = as.numeric(time_day))   # numeric for continuous color scale
arrow_step <- 100   # every Nth point will become an arrow; decrease for denser arrows
arrow_pts <- df_plot %>%
  slice(seq(1, n(), by = arrow_step)) %>%
  mutate(xend = lead(Shallow), yend = lead(Deep)) %>%
  drop_na()

time_breaks_num <- as.numeric(c(min(df_plot$time_num),
                                median(df_plot$time_num),
                                max(df_plot$time_num)))
time_breaks_label <- as_datetime(time_breaks_num) %>% format("%Y")

p5 <- ggplot(df_plot, aes(Shallow, Deep, color = time_num)) +
  geom_path(alpha = 0.45, size = 0.35) +
  geom_segment(data = arrow_pts,
               aes(x = Shallow, y = Deep, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
               inherit.aes = FALSE,
               color = "grey20", size = 0.35, alpha = 0.9) +
  # 1:1 reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  scale_color_viridis_c(option = "plasma",
                        name = "Date",
                        breaks = time_breaks_num,
                        labels = time_breaks_label) +
  coord_equal(expand = TRUE) +
  labs(
    x = "Shallow Temperature (°C)",
    y = "Deep Temperature (°C)",
  ) +
  theme(
    legend.position = c(0.92, 0.5), 
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(1, "cm") 
  ) +
  theme_pubr(base_size = 12)
p5

####################### Density plot by month

df_month_ridges <- df_long_10min %>%
  mutate(month = month(time_10min, label = TRUE)) 

p6 <- ggplot(df_month_ridges, aes(Temp, month, fill = Depth)) +
  geom_density_ridges(alpha = 0.6) +
  theme_pubr() + scale_fill_manual(values = c("Deep" = "#00AFBB", "Shallow" = "#E7B800")) +
  theme_cowplot() +
  labs(x = "Temperature (°C)", 
       y = "Month") +
  theme(legend.position = "none")
p6
# 
# ## plot with box
# p6 <- ggplot(df_month_ridges, aes(Temp, month, fill = Depth)) +
#   geom_density_ridges(alpha = 0.6) +
#   annotate("rect",
#            xmin = 28,
#            xmax = 33,
#            ymin = 4.8,   # Just below June (6)
#            ymax = 9.4,   # Just above August (8)
#            color = "steelblue",
#            fill = NA,    # Keep it transparent inside
#            linewidth = 2) +
#   theme_cowplot() +
#   scale_fill_manual(values = c("Deep" = "#00AFBB", "Shallow" = "#E7B800")) +
#   labs(x = "Temperature (°C)", y = "Month") +
#   theme(legend.position = "none")
# p6

###########################

## combined
p2_leg <- p2 + 
  theme(
    legend.position = "right",
    legend.background = element_rect(
      color = "black",      
      fill = "white",       
      linewidth = 0.5,      
      linetype = "solid"    
    ),
    legend.margin = margin(3, 3, 3, 3) # Removed the + from inside here
  ) + 
  labs(color = "Site Depth") 

shared_legend <- get_legend(p2_leg)
p1_noleg <- p1 + theme(legend.position = "none")
p1_noleg <- p1_noleg + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = -10))
p2_noleg <- p2 + theme(legend.position = "none")
p3_noleg <- p3 + theme(legend.position = "none")

top_row <- plot_grid(
  p2_noleg, p4, p3_noleg,
  labels = c("a", "b", "c"),
  label_size = 20,
  label_fontface = "bold",
  rel_widths = c(1.1, 1.1, 0.7),
  nrow = 1,
  align = "hv"
)
bottom_plots <- plot_grid(
  p1_noleg, p6, p5, 
  labels = c("d", "e", "f"),
  label_size = 20,
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1, 1.3, 0.8), 
  align = "v",     
  axis = "t"
)
bottom_row <- plot_grid(
  shared_legend, 
  bottom_plots, 
  nrow = 1, 
  rel_widths = c(0.12, 0.88)
)

# Combine everything + legend
final_plot <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1)
)
ggsave("figs/combined_alldata.jpg", final_plot, width=12, height=7, unit="in")

################################################################################

##hourly autocorrelation analysis
shallow_hr <- subset(df_long_hour, Depth == "Shallow")
deep_hr <- subset(df_long_hour, Depth == "Deep")

shallow_hr.auto <- acf(shallow_hr$Temp, plot = FALSE)
deep_hr.auto <- acf(deep_hr$Temp, plot = FALSE)

shallow_hr.auto <- as.data.frame(shallow_hr.auto$acf[0:47])
shallow_hr.auto$Lag <- seq(from = 1, to = 47, by = 1)
colnames(shallow_hr.auto) <- c("Autocorr", "Lag")

deep_hr.auto <- as.data.frame(deep_hr.auto$acf[0:47])
deep_hr.auto$Lag <- seq(from = 1, to = 47, by = 1)
colnames(deep_hr.auto) <- c("Autocorr", "Lag")

shallow_hr.auto$Depth <- "Shallow"
deep_hr.auto$Depth <- "Deep"

all.auto <- rbind(shallow_hr.auto, deep_hr.auto)

#plotting autocorrelation with hourly lag
p21 <- ggplot(all.auto, aes(x = Lag, y = Autocorr, shape = Depth, colour = Depth)) +
  geom_point(size = 2)+
  scale_shape_manual("Site Depth", values = c("Shallow" = 16, "Deep" = 18),
                     breaks = c("Shallow", "Deep")) +
  scale_colour_manual("Site Depth", values = c("Shallow" = "#E7B800","Deep" = "#00AFBB"),
                      breaks = c("Shallow", "Deep")) +
  ylab("Temperature Autocorrelation") +
  xlab("Lag (hours)") +
  theme_pubr() + 
  theme(axis.title = element_text(size = 16))
p21

##daily autocorrelation analysis
shallow_day <- subset(df_long_day, Depth == "Shallow")
deep_day <- subset(df_long_day, Depth == "Deep")

shallow_day.auto <- acf(shallow_day$Temp, plot = FALSE)
deep_day.auto <- acf(deep_day$Temp, plot = FALSE)

shallow_day.auto <- as.data.frame(shallow_day.auto$acf[0:33])
shallow_day.auto$Lag <- seq(from = 1, to = 33, by = 1)
colnames(shallow_day.auto) <- c("Autocorr", "Lag")

deep_day.auto <- as.data.frame(deep_day.auto$acf[0:33])
deep_day.auto$Lag <- seq(from = 1, to = 33, by = 1)
colnames(deep_day.auto) <- c("Autocorr", "Lag")

shallow_day.auto$Depth <- "Shallow"
deep_day.auto$Depth <- "Deep"

all_day.auto <- rbind(shallow_day.auto, deep_day.auto)

#plotting autocorrelation with daily lag
p22 <- ggplot(all_day.auto, aes(x = Lag, y = Autocorr, shape = Depth, colour = Depth)) +
  geom_point(size = 2)+
  scale_shape_manual("Site Depth", values = c("Shallow" = 16, "Deep" = 18),
                     breaks = c("Shallow", "Deep")) +
  scale_colour_manual("Site Depth", values = c("Shallow" = "#E7B800","Deep" = "#00AFBB"),
                      breaks = c("Shallow", "Deep")) +
  ylab("Temperature Autocorrelation") +
  xlab("Lag (days)") +
  theme_pubr() +
  theme(axis.title = element_text(size = 16))
p22


##combine into singular plot
p21_noleg <- p21 + theme(legend.position = "none")
p22_leg <- p22 + theme(
  legend.position = "top",
  legend.direction = "horizontal", 
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16),
  legend.background = element_blank(), 
  legend.key = element_blank(),     
  legend.margin = margin(b = -10)     
)

shared_legend_2 <- get_legend(p22_leg)
p21_noleg <- p21_noleg + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = -10))
p22_noleg <- p22_leg + theme(legend.position="none")

top_row_2 <- plot_grid(
  shared_legend_2, 
  nrow = 1
)
bot_row_2 <- plot_grid(
  p21_noleg, p22_noleg,
  labels = c("a", "b"),
  label_size = 20,
  label_fontface = "bold",
  rel_widths = c(1,1),
  nrow = 1,
  align = "hv"
)

final_plot_2 <- plot_grid(
  top_row_2,
  bot_row_2,
  ncol = 1,
  rel_heights = c(0.08, 1)
)
ggsave("figs/acf_all_data.jpg", final_plot_2, width=12, height=7, unit="in")

