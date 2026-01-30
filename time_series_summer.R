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
library(scales)
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

df_long <- df_long %>%
  filter(month(time_clean) %in% c(5, 6, 7, 8))

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

mean(shallow$Temp) # 29.89178
range(shallow$Temp) # 25.65707 - 33.99100
sd(shallow$Temp) # 1.390386

mean(deep$Temp) # 29.61793
range(deep$Temp) # 25.16131 - 32.50217
sd(deep$Temp) # 1.408201

#KS test
ks.test(shallow$Temp,deep$Temp) 
# D = 0.068945, p-value < 2.2e-16

#Permanova on temperature distributions
perm.anova(Temp ~ Depth, data = df_long_10min, nperm = 999)
# 999 permutations
# Sum Sq    Df Mean Sq F value Pr(>F)    
# Depth       1832     1 1832.28  935.75  0.001 ***
#   Residuals 191361 97728    1.96                   

#Check equality of variances
leveneTest(Temp ~ Depth, data = df_long_10min, na.action = na.pass)
# Levene's Test for Homogeneity of Variance (center = median: na.pass)
#          Df F value    Pr(>F)    
# group     1  41.837 9.966e-11 ***
#       97728                                       

#Percentage of time above 30C
# https://floridakeys.noaa.gov/coral-bleaching-faqs.html
df_long_10min %>%
  dplyr::mutate(exceed = Temp > 30) %>%
  dplyr::group_by(Depth) %>%
  dplyr::summarize(pct_time = mean(exceed) * 100)
# Deep 48
# Shallow 53.7

################################################################################

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
df_faceted <- df_long_10min %>%
  mutate(year_label = year(time_10min))

p2 <- ggplot(df_faceted, aes(x = time_10min, y = Temp, color = Depth)) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("Deep" = "#00AFBB", "Shallow" = "#E7B800")) +
  theme_cowplot() +
  # Stack vertically (1 column) and allow x-axes to vary by year
  facet_wrap(~year_label, scales = "free_x", ncol = 1) +
  # Customize X-axis ticks
  scale_x_datetime(
    breaks = date_breaks("1 month"), # Ensure a break exists for every month
    labels = label_date("%b")) + 
  labs(x = "Month", y = "Temperature (°C)") +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 13), 
    axis.text = element_text(size = 11),
    # Rotate x-axis ticks by 48 degrees
    axis.text.x = element_text(angle = 48, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

final_p2 <- plot_grid(p2) 
final_faceted_plot <- ggdraw(final_p2) +
  draw_label(
    "All KS p < 0.001", 
    x = 0.02,           
    y = 0.01,          
    hjust = 0,         
    vjust = 0,        
    size = 12, 
    fontface = "italic"
  )
final_faceted_plot
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
  scale_x_continuous(limits = c(28, 34), breaks = c(28, 30, 32, 34)) +
  scale_y_continuous(limits = c(28, 34), breaks = c(28, 30, 32, 34)) +
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
# D = 0.57018, p-value < 2.2e-16

########## Box plots
df_june <- df_long_day %>% filter(month(time_day) == 6)
df_july <- df_long_day %>% filter(month(time_day) == 7)
df_august <- df_long_day %>% filter(month(time_day) == 8)

# t-test p = 0.32
p6 <- ggboxplot(df_august, x = "Depth", y = "Temp",
             fill = "Depth", palette = c("#E7B800","#00AFBB"),
             add = "jitter", 
             add.params = list(alpha = 0.5, size = 1.5)) + 
  scale_y_continuous(limits = c(28, 33.5)) +
  annotate("text",
           x = -Inf, y = Inf,             # Top Left corner
           label = "M-W U-test, p = 0.32",
           hjust = -0.1, vjust = 1.5,      # Nudge it inside the plot
           size = 3.5,
           fontface = "plain") +
  labs(y = "August Temperature (°C)", x = "") +
  theme_cowplot() +
  theme(legend.position = "none")
p6

# p6 <- p6 +
#   stat_compare_means(
#     method = "wilcox.test",
#     label.x.npc = "left",
#     label.y.npc = "top",
#     aes(label = paste0("Mann-Whitney U, ", ..p.format..)),
#     size = 3.5
#   )

## t-test p = 0.025
p4 <- ggboxplot(df_july, x = "Depth", y = "Temp",
                fill = "Depth", palette = c("#E7B800","#00AFBB"),
                add = "jitter", 
                add.params = list(alpha = 0.5, size = 1.5)) + 
  scale_y_continuous(limits = c(28, 33.5)) +
  annotate("text",
           x = -Inf, y = Inf,             # Top Left corner
           label = "M-W U-test, p = 0.0066",
           hjust = -0.1, vjust = 1.5,      # Nudge it inside the plot
           size = 3.5,
           fontface = "plain") +
  labs(y = "July Temperature (°C)", x = "") +
  theme_cowplot() +
  theme(legend.position = "none")
p4
 
# p4 <- p4 + 
#   stat_compare_means(
#     method = "wilcox.test", 
#     label.x.npc = "left", 
#     label.y.npc = "top",
#     aes(label = paste0("Mann-Whitney U, ", ..p.format..)),
#     size = 3.5
#   )

## t-test p = 0.0079
p5 <- ggboxplot(df_june, x = "Depth", y = "Temp",
                fill = "Depth", palette = c("#E7B800","#00AFBB"),
                add = "jitter", add.params = list(alpha = 0.4)) +
  scale_y_continuous(limits = c(28, 33.5)) +
  annotate("text",
           x = -Inf, y = Inf,             # Top Left corner
           label = "M-W U-test, p = 0.00076",
           hjust = -0.1, vjust = 1.5,      # Nudge it inside the plot
           size = 3.5,
           fontface = "plain") +
  labs(y = "June Temperature (°C)", x = "") +
  theme_cowplot() +
  theme(legend.position = "none")
p5 

# p5 <- p5 +
#   stat_compare_means(
#     method = "wilcox.test",
#     label.x.npc = "left",
#     label.y.npc = "top",
#     aes(label = paste0("Mann-Whitney U, ", ..p.format..)),
#     size = 3.5
#   )
# p5

################## Combine figures

p2_leg <- p2 + theme(
  legend.position = "top",
  legend.direction = "horizontal", 
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16),
  legend.key = element_blank(),     
  # legend.margin = margin(b = -10),
  legend.background = element_rect(
      color = "black",      
      fill = "white",       
      linewidth = 0.5,      
      linetype = "solid"    
    ),
  legend.margin = margin(3, 3, 3, 3)
) + labs(color = "Site Depth") 

shared_legend <- get_legend(p2_leg)
p1_noleg <- p1 + theme(legend.position = "none")

left_top_row <- plot_grid(
  p1_noleg, p3, 
  labels = c("a", "b"), 
  label_size = 20, 
  label_fontface = "bold",
  nrow = 1,
  align = "h"
)

legend_row <- plot_grid(
  NULL, shared_legend, NULL, 
  ncol = 3, 
  rel_widths = c(1.5, 2, 1) # Adjust these values to shift the legend left or right
)

left_bottom_row <- plot_grid(
  p5, p4, p6, 
  labels = c("c", "d", "e"), 
  label_size = 20, 
  label_fontface = "bold",
  nrow = 1,
  align = "h",
  axis = "bt" # Align bottoms and tops of the boxplot axes
)

left_column <- plot_grid(
  legend_row, 
  left_top_row, 
  left_bottom_row, 
  ncol = 1, 
  rel_heights = c(0.15, 1, 1) 
)

final_big_plot <- plot_grid(
  left_column, 
  final_faceted_plot, 
  labels = c("", "f"), 
  label_size = 20,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1.5, 1)
)
ggsave("figs/summer_data.jpg", final_big_plot, width=18, height=10, unit="in")


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
ggsave("figs/acf_summerdata.jpg", final_plot_2, width=12, height=7, unit="in")