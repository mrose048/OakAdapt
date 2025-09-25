##%######################################################%##
#                                                          #
####            Dye Creek oak change results            ####
#                                                          #
##%######################################################%##

library(caTools)
library(reticulate)
require(terra)
require(dplyr)
require(sf)
require(elevatr)
require(tidyverse)
require(ggspatial)
library(tidyr)
library(shiny)

# setwd('H:/My Drive/16_Oak_Woodland/1_Data/')

dyecreek <- vect('1_Data/00_Shapefiles/study_areas/DyeCreek.shp') 


##%######################################################%##
#                                                          #
####                      Landsat                       ####
#                                                          #
##%######################################################%##

# Set file paths (adjust as needed)
veg_path <- "1_Data/Landsat/05_DyeCreek_temporal_stack/DyeCreek_LS_1982_2025_veg_stack"
tmr_pca <- "1_Data/Landsat/06_DyeCreek_tem_unmixed/full/DyeCreek_LS_1982_2025_TMR_PCA"

# Load raster stacks
veg_stack <- rast(veg_path)
tmr_pca <- rast(tmr_pca) 
tmr_pca[[1:3]]

# Create sample points within study area
set.seed(2024)
sample_points <- spatSample(veg_stack[[1]], size = 1000, method = "random", as.points = TRUE)

# Extract veg time series and TMR values at points
veg_vals <- terra::extract(veg_stack, sample_points)
tmr_vals <- terra::extract(tmr_pca, sample_points)

# Combine into data frame
veg_tmr_sample_df <- cbind(veg_vals, tmr_vals[, 2:4]) # remove ID column in tmr_vals

# Rename and clean
colnames(veg_tmr_sample_df)[(ncol(veg_tmr_sample_df)-2):ncol(veg_tmr_sample_df)] <- c("TMR1", "TMR2", "TMR3")
veg_tmr_sample_df$ID <- NULL

# Add pixel ID
veg_tmr_sample_df$pixel_id <- seq_len(nrow(veg_tmr_sample_df))

# Pivot to long format
veg_tmr_sample_df_long <- veg_tmr_sample_df |>
  pivot_longer(
    cols = matches("^[0-9]{4}"),  # columns that start with a 4-digit year
    names_to = "date",
    values_to = "veg_score"
  ) |>
  mutate(
    date = as.numeric(date),
    TMR1_cat = case_when(
      TMR1 >= quantile(TMR1, 0.95, na.rm = TRUE) ~ "TMR1_high",
      TMR1 <= quantile(TMR1, 0.05, na.rm = TRUE) ~ "TMR1_low",
      TRUE ~ "TMR1_mid"
    ),
    TMR2_cat = case_when(
      TMR2 >= quantile(TMR2, 0.9, na.rm = TRUE) ~ "TMR2_high",
      TMR2 <= quantile(TMR2, 0.1, na.rm = TRUE) ~ "TMR2_low",
      TRUE ~ "TMR2_mid"
    ),
    TMR3_cat = case_when(
      TMR3 >= quantile(TMR3, 0.95, na.rm = TRUE) ~ "TMR3_high",
      TMR3 <= quantile(TMR3, 0.05, na.rm = TRUE) ~ "TMR3_low",
      TRUE ~ "TMR3_mid"
    )
  )


# Compute mean and sd of veg_score by year and TMR category
avg_trends <- bind_rows(
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR1_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR2_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR3_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE))
)


ui <- fluidPage(
  titlePanel("Vegetation Trends by TMR Category"),
  sidebarLayout(
    sidebarPanel(
      selectInput("tmr_component", "Choose TMR Component:",
                  choices = c("TMR1", "TMR2", "TMR3"),
                  selected = "TMR3"),
      checkboxGroupInput("tmr_level", "TMR Level:",
                         choices = c("High" = "high", "Low" = "low"),
                         selected = c("high", "low")),
      sliderInput("y_range", "Y-axis Range:", min = -0.2, max = 1, value = c(-0.1, 0.6), step = 0.05),
      checkboxInput("add_vlines", "Show Event Years", value = TRUE)
    ),
    mainPanel(
      plotOutput("trendPlot", height = "600px")
    )
  )
)

server <- function(input, output) {
  output$trendPlot <- renderPlot({
    selected_categories <- paste0(input$tmr_component, "_", input$tmr_level)

    plot_data <- avg_trends |>
      filter(TMR_category %in% selected_categories)

    p <- ggplot(plot_data, aes(x = date, y = mean_v_score, color = TMR_category)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(
        "TMR1_high" = "red", "TMR1_low" = "blue",
        "TMR2_high" = "darkgreen", "TMR2_low" = "lightgreen",
        "TMR3_high" = "purple", "TMR3_low" = "orange"
      )) +
      theme_minimal(base_size = 16) +
      labs(
        title = paste("Vegetation Trends by", input$tmr_component),
        x = "Year", y = "Mean Vegetation Fraction",
        color = "Category"
      ) +
      ylim(input$y_range)

    if (input$add_vlines) {
      event_years <- c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025)
      p <- p + geom_vline(xintercept = event_years, linetype = "dashed", color = "gray70")
    }

    p
  })
}

shinyApp(ui = ui, server = server)

veg_summary <- avg_trends %>%
  filter(date <= 2024.45) %>%
  mutate(period = ifelse(date <= 1987, "Early", "Late")) %>%
  group_by(TMR_category, period) %>%
  summarise(mean_v = mean(mean_v_score, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = mean_v) %>%
  mutate(
    prop_change = (Late - Early) / Early,
    pct_change = 100 * prop_change
  )

# TMR 3 trend
tmr3_trend <- veg_tmr_sample_df_long |>
  mutate(
    TMR3_group = case_when(
      TMR3 <= quantile(TMR3, 0.05, na.rm = TRUE) ~ "Q1 (lowest)",
      TMR3 > quantile(TMR3, 0.05, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.35, na.rm = TRUE) ~ "Q2",
      TMR3 > quantile(TMR3, 0.35, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.65, na.rm = TRUE) ~ "Q3",
      TMR3 > quantile(TMR3, 0.65, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.95, na.rm = TRUE) ~ "Q4",
      TMR3 > quantile(TMR3, 0.95, na.rm = TRUE) ~ "Q5 (highest)",
      TRUE ~ NA_character_
    ),
    TMR3_group = factor(TMR3_group, levels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
  )

# TMR 1 trend
tmr1_trend <- veg_tmr_sample_df_long |>
  mutate(
    TMR1_group = case_when(
      TMR1 <= quantile(TMR1, 0.05, na.rm = TRUE) ~ "Q1 (lowest)",
      TMR1 > quantile(TMR1, 0.1, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.35, na.rm = TRUE) ~ "Q2",
      TMR1 > quantile(TMR1, 0.35, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.65, na.rm = TRUE) ~ "Q3",
      TMR1 > quantile(TMR1, 0.65, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.9, na.rm = TRUE) ~ "Q4",
      TMR1 > quantile(TMR1, 0.9, na.rm = TRUE) ~ "Q5 (highest)",
      TRUE ~ NA_character_
    ),
    TMR1_group = factor(TMR1_group, levels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
  )

quintile_trends <- tmr3_trend |>
  group_by(date, TMR3_group) |>
  summarise(
    mean_v_score = mean(veg_score, na.rm = TRUE),
    .groups = "drop"
  )


p_quintiles3 <- ggplot(
  quintile_trends %>% 
    na.omit() %>% 
    filter(TMR3_group %in% c('Q1 (lowest)', 'Q5 (highest)')),
  aes(x = date, y = mean_v_score, color = TMR3_group)
) +
  geom_line(size = 0.5, alpha = 0.9) +  # raw fluctuations
  geom_smooth(se = FALSE, size = 1.5) + # smoothed trend
  scale_color_manual(values = c(
    "Q1 (lowest)" = "#0D0887FF",
    "Q5 (highest)" = "#B8860B"
  )) +
  labs(
    title = "",
    x = "Year",
    y = "Mean Vegetation Fraction",
    color = "TMR3 Quintile"
  ) +
  theme_minimal(base_size = 26) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26),
    legend.position = 'bottom'
  ) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  annotate("rect", xmin = 1987, xmax = 1992, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
  annotate("rect", xmin = 2012, xmax = 2015.7, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
  annotate("rect", xmin = 2020, xmax = 2022, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
  annotate("rect", xmin = 2007, xmax = 2009, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15)
print(p_quintiles3)

# save
ggsave(
  "3_Figures/Manuscript_figures/DyeCreek_LS_TMR3_residuals_quintiles.png",
  plot = p_quintiles3,
  width = 15,
  height = 10,
  dpi = 600
)

quintile_trends <- tmr1_trend |>
  group_by(date, TMR1_group) |>
  summarise(
    mean_v_score = mean(veg_score, na.rm = TRUE),
    .groups = "drop"
  )

p_quintiles <- ggplot(
  quintile_trends %>% na.omit() %>% filter(date >= 2010) %>% filter(TMR1_group %in% c('Q1 (lowest)', 'Q5 (highest)')),
  aes(x = date, y = mean_v_score, color = TMR1_group)
) +
  geom_line(size = 1) +
  scale_color_manual(values = c(
    "Q1 (lowest)" = "#0D0887FF",
    # "Q2" = "#7301A8FF",
    # "Q3" = "#BD3786FF",
    # "Q4" = "#ED7953FF",
    "Q5 (highest)" = "#FDC926FF"
  )) +
  labs(
    title = "Vegetation Trends by TMR1 Quintile",
    x = "Year",
    y = "Mean Vegetation Fraction",
    color = "TMR3 Quintile"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    legend.position = 'none',
    legend.background = element_rect(fill = "white", color = "black"),
    legend.direction = 'horizontal'
  ) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025)) +
  # annotate(
  #   "rect",
  #   xmin = 1987,
  #   xmax = 1992,
  #   ymin = -Inf,
  #   ymax = Inf,
  #   fill = "red",
  #   alpha = 0.25
  # ) +
  annotate(
    "rect",
    xmin = 2012,
    xmax = 2015.7,
    ymin = -Inf,
    ymax = Inf,
    fill = "red",
    alpha = 0.25
  ) +
  annotate(
    "rect",
    xmin = 2020,
    xmax = 2022,
    ymin = -Inf,
    ymax = Inf,
    fill = "red",
    alpha = 0.25
  )

print(p_quintiles)

# save
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/DyeCreek_LS_TMR1_residuals_quintiles.png",
  plot = p_quintiles,
  width = 15,
  height = 12,
  dpi = 600
)

# Prepare raster
r <- tmr_pca[[3]]
raster_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
colnames(raster_df) <- c("x", "y", "value")

# Define uneven breaks once, based on TMR3 values (raster or sample points)
uneven_breaks <- quantile(tmr_pca[[3]][], probs = c(0, 0.05, 0.35, 0.65, 0.95, 1), na.rm = TRUE)
uneven_labels <- c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")

# For raster data
raster_df <- as.data.frame(tmr_pca[[3]], xy = TRUE, na.rm = TRUE)
colnames(raster_df) <- c("x", "y", "value")

raster_df <- raster_df %>%
  mutate(
    quintile = cut(value,
                   breaks = uneven_breaks,
                   include.lowest = TRUE,
                   labels = uneven_labels)
  )


# Assign reversed plasma colors
plasma_colors <- c("#0D0887FF","#7301A8FF", "#BD3786FF",  "#ED7953FF","#FDC926FF")
names(plasma_colors) <- uneven_labels

legend_title <- "PC3 Quintile"
bbox <- st_bbox(dyecreek)

# Plot
p <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = quintile)) +
  scale_fill_manual(values = plasma_colors, name = legend_title) +
  geom_sf(data = st_as_sf(dyecreek), fill = NA, color = "black", lwd = 1) +
  annotation_scale(location = "bl", width_hint = 0.17, height = unit(0.6, "cm"), text_cex = 1.5) +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(2.5, "cm"), width = unit(1.75, "cm")) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.position = c(.8, .1),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )

p


library(patchwork)

# Combine the plots side-by-side or top-bottom
combined_plot <- (p_quintiles + p) +
  plot_layout(widths = c(0.4, 0.6)) 

print(combined_plot)

# save
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/DyeCreek_LS_TMR3_residuals_quintiles.png",
  plot = combined_plot,
  width = 30,
  height = 12,
  dpi = 600
)

# save 
writeRaster(tmr_pca[[3]], "08_Outputs/DyeCreek_Landsat_change_raw.tif") #high values = resilience, low values = struggling
writeRaster(tmr_pca[[1]], "08_Outputs/DyeCreek_Landsat_change_raw_TMR1.tif") #high values = resilience, low values = struggling

##%######################################################%##
#                                                          #
####                     Sentinel-2                     ####
#                                                          #
##%######################################################%##

# Set file paths (adjust as needed) 
veg_path <- "1_Data/Sentinel_2/05_DyeCreek_temporal_stack/DyeCreek_S2_stack"
# tmr_pca <- "Sentinel2_folder/06_DyeCreek_tem_unmixed/full/DyeCreek_S2_TMR_full_PCA"
tmr_pca <- "1_Data/Sentinel_2/06_DyeCreek_tem_unmixed/full/DyeCreek_S2_TMR_full_PCA"

# Load raster stacks
veg_stack <- rast(veg_path)
tmr_pca <- rast(tmr_pca) 
tmr_pca[[1:3]]

# Extract the third principal component layer
tmr3 <- tmr_pca[[3]]

# Create a mask of cells that are not NA
non_na_cells <- which(!is.na(values(tmr3)))

# Use those indices to sample from only valid cells
# Number of samples you want:
n_samples <- 2000

# Randomly sample indices from valid cells
set.seed(123)  # for reproducibility
sample_indices <- sample(non_na_cells, size = n_samples)

# Convert cell indices to xy coordinates
sample_points <- xyFromCell(tmr3, sample_indices)

# Convert to SpatVector
sample_points_vect <- vect(sample_points, type = "points", crs = crs(tmr3))


# Extract veg time series and TMR values at points
veg_vals <- terra::extract(veg_stack, sample_points_vect)
tmr_vals <- terra::extract(tmr_pca, sample_points_vect)

# Combine into data frame
veg_tmr_sample_df <- cbind(veg_vals, tmr_vals[, 2:4]) # remove ID column in tmr_vals

# Rename and clean
colnames(veg_tmr_sample_df)[(ncol(veg_tmr_sample_df)-2):ncol(veg_tmr_sample_df)] <- c("TMR1", "TMR2", "TMR3")
veg_tmr_sample_df$ID <- NULL

# Add pixel ID
veg_tmr_sample_df$pixel_id <- seq_len(nrow(veg_tmr_sample_df))

# Pivot to long format
veg_tmr_sample_df_long <- veg_tmr_sample_df |>
  pivot_longer(
    cols = matches("^[0-9]{4}"),  # columns that start with a 4-digit year
    names_to = "date",
    values_to = "veg_score"
  ) |>
  mutate(
    date = as.numeric(date),
    TMR1_cat = case_when(
      TMR1 >= quantile(TMR1, 0.95, na.rm = TRUE) ~ "TMR1_high",
      TMR1 <= quantile(TMR1, 0.05, na.rm = TRUE) ~ "TMR1_low",
      TRUE ~ "TMR1_mid"
    ),
    TMR2_cat = case_when(
      TMR2 >= quantile(TMR2, 0.95, na.rm = TRUE) ~ "TMR2_high",
      TMR2 <= quantile(TMR2, 0.05, na.rm = TRUE) ~ "TMR2_low",
      TRUE ~ "TMR2_mid"
    ),
    TMR3_cat = case_when(
      TMR3 >= quantile(TMR3, 0.95, na.rm = TRUE) ~ "TMR3_high",
      TMR3 <= quantile(TMR3, 0.05, na.rm = TRUE) ~ "TMR3_low",
      TRUE ~ "TMR3_mid"
    )
  )


# Compute mean and sd of veg_score by year and TMR category
avg_trends <- bind_rows(
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR1_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR2_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE)),
  veg_tmr_sample_df_long |> group_by(date, TMR_category = TMR3_cat) |> summarise(mean_v_score = mean(veg_score, na.rm = TRUE))
)


ui <- fluidPage(
  titlePanel("Vegetation Trends by TMR Category"),
  sidebarLayout(
    sidebarPanel(
      selectInput("tmr_component", "Choose TMR Component:",
                  choices = c("TMR1", "TMR2", "TMR3"),
                  selected = "TMR3"),
      checkboxGroupInput("tmr_level", "TMR Level:",
                         choices = c("High" = "high", "Low" = "low"),
                         selected = c("high", "low")),
      sliderInput("y_range", "Y-axis Range:", min = -0.2, max = 1, value = c(-0.1, 0.6), step = 0.05),
      checkboxInput("add_vlines", "Show Event Years", value = TRUE)
    ),
    mainPanel(
      plotOutput("trendPlot", height = "600px")
    )
  )
)

server <- function(input, output) {
  output$trendPlot <- renderPlot({
    selected_categories <- paste0(input$tmr_component, "_", input$tmr_level)

    plot_data <- avg_trends |>
      filter(TMR_category %in% selected_categories)

    p <- ggplot(plot_data, aes(x = date, y = mean_v_score, color = TMR_category)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(
        "TMR1_high" = "red", "TMR1_low" = "blue",
        "TMR2_high" = "darkgreen", "TMR2_low" = "lightgreen",
        "TMR3_high" = "purple", "TMR3_low" = "orange"
      )) +
      theme_minimal(base_size = 16) +
      labs(
        title = paste("Vegetation Trends by", input$tmr_component),
        x = "Year", y = "Mean Vegetation Fraction",
        color = "Category"
      ) +
      ylim(input$y_range)

    if (input$add_vlines) {
      event_years <- c(2016, 2018, 2020, 2022, 2024)
      p <- p + geom_vline(xintercept = event_years, linetype = "dashed", color = "gray70")
    }

    p
  })
}

shinyApp(ui = ui, server = server)

# TMR 3 trend
tmr3_trend <- veg_tmr_sample_df_long |>
  mutate(
    TMR3_group = case_when(
      TMR3 <= quantile(TMR3, 0.05, na.rm = TRUE) ~ "Q1 (lowest)",
      TMR3 > quantile(TMR3, 0.05, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.35, na.rm = TRUE) ~ "Q2",
      TMR3 > quantile(TMR3, 0.35, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.65, na.rm = TRUE) ~ "Q3",
      TMR3 > quantile(TMR3, 0.65, na.rm = TRUE) & TMR3 <= quantile(TMR3, 0.95, na.rm = TRUE) ~ "Q4",
      TMR3 > quantile(TMR3, 0.95, na.rm = TRUE) ~ "Q5 (highest)",
      TRUE ~ NA_character_
    ),
    TMR3_group = factor(TMR3_group, levels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
  )

# TMR 1 trend
tmr1_trend <- veg_tmr_sample_df_long |>
  mutate(
    TMR1_group = case_when(
      TMR1 <= quantile(TMR1, 0.05, na.rm = TRUE) ~ "Q1 (lowest)",
      TMR1 > quantile(TMR1, 0.05, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.35, na.rm = TRUE) ~ "Q2",
      TMR1 > quantile(TMR1, 0.35, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.65, na.rm = TRUE) ~ "Q3",
      TMR1 > quantile(TMR1, 0.65, na.rm = TRUE) & TMR1 <= quantile(TMR1, 0.95, na.rm = TRUE) ~ "Q4",
      TMR1 > quantile(TMR1, 0.95, na.rm = TRUE) ~ "Q5 (highest)",
      TRUE ~ NA_character_
    ),
    TMR1_group = factor(TMR1_group, levels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
  )


quintile_trends <- tmr3_trend |>
  group_by(date, TMR3_group) |>
  summarise(
    mean_v_score = mean(veg_score, na.rm = TRUE),
    .groups = "drop"
  )

# quintile_trends <- tmr1_trend |>
#   group_by(date, TMR1_group) |>
#   summarise(
#     mean_v_score = mean(veg_score, na.rm = TRUE),
#     .groups = "drop"
#   )

p_quintiles3 <- ggplot(
  quintile_trends %>%
    na.omit() %>%
    filter(TMR3_group %in% c('Q1 (lowest)', 'Q5 (highest)')),
  aes(x = date, y = mean_v_score, color = TMR3_group)
) +
  geom_line(size = 0.5, alpha = 0.9) +  # raw fluctuations
  geom_smooth(se = FALSE, size = 1.5) + # smoothed trend
  scale_color_manual(values = c(
    "Q1 (lowest)" = "#0D0887FF",
    "Q5 (highest)" = "#B8860B"
  )) +
  labs(
    title = "",
    x = "Year",
    y = "Mean Vegetation Fraction",
    color = "TMR3 Quintile"
  ) +
  theme_minimal(base_size = 26) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26),
    legend.position = 'bottom'
  ) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  annotate("rect", xmin = 1984, xmax = 1985, ymin = -Inf, ymax = Inf, fill = "transparent", alpha = 0.15) +
  annotate("rect", xmin = 1987, xmax = 1992, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
  annotate("rect", xmin = 2012, xmax = 2015.7, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
  annotate("rect", xmin = 2020, xmax = 2022, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
  annotate("rect", xmin = 2007, xmax = 2009, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15)
print(p_quintiles3)

# save
ggsave(
  "3_Figures/Manuscript_figures/DyeCreek_S2_TMR3_residuals_quintiles.png",
  plot = p_quintiles3,
  width = 15,
  height = 10,
  dpi = 600
)


# Prepare raster
r <- tmr_pca[[2]]
raster_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
colnames(raster_df) <- c("x", "y", "value")

# Define uneven breaks once, based on TMR3 values (raster or sample points)
uneven_breaks <- quantile(tmr_pca[[2]][], probs = c(0, 0.05, 0.35, 0.65, 0.95, 1), na.rm = TRUE)
uneven_labels <- c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")

# For raster data
raster_df <- as.data.frame(tmr_pca[[2]], xy = TRUE, na.rm = TRUE)
colnames(raster_df) <- c("x", "y", "value")

raster_df <- raster_df %>%
  mutate(
    quintile = cut(value,
                   breaks = uneven_breaks,
                   include.lowest = TRUE,
                   labels = uneven_labels)
  )


# Assign reversed plasma colors
plasma_colors <- c("#0D0887FF","#7301A8FF", "#BD3786FF",  "#ED7953FF","#FDC926FF")
names(plasma_colors) <- uneven_labels

legend_title <- "PC2 Quintile"
bbox <- st_bbox(dyecreek)

# Plot
p <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = quintile)) +
  scale_fill_manual(values = plasma_colors, name = legend_title) +
  geom_sf(data = st_as_sf(dyecreek), fill = NA, color = "black", lwd = 1) +
  annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.6, "cm"), text_cex = 1.5) +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(2.5, "cm"), width = unit(1.75, "cm")) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.position = c(.8, .1),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )

p


library(patchwork)

# Combine the plots side-by-side or top-bottom
combined_plot <- (p_quintiles + p) +
  plot_layout(widths = c(0.4, 0.6)) 

print(combined_plot)

# save
ggsave(
  "H:/My Drive/16_Oak_Woodland/3_Figures/Manuscript_figures/RandallLS_TMR3_residuals_quintiles_pre_fire.png",
  plot = combined_plot,
  width = 30,
  height = 12,
  dpi = 600
)


library(terra)

# Select the raster layer
r <- tmr_pca[[3]] 
r1 <- tmr_pca[[1]] *-1

# Define custom quantile breaks
uneven_breaks <- quantile(r[], probs = c(0, 0.05, 0.35, 0.65, 0.95, 1), na.rm = TRUE)

# Reclassify raster based on breaks into categories 1–5
r_class <- classify(r, rcl = matrix(c(
  uneven_breaks[1], uneven_breaks[2], 1,
  uneven_breaks[2], uneven_breaks[3], 2,
  uneven_breaks[3], uneven_breaks[4], 3,
  uneven_breaks[4], uneven_breaks[5], 4,
  uneven_breaks[5], uneven_breaks[6], 5
), ncol = 3, byrow = TRUE))

# Assign factor labels
levels(r_class) <- data.frame(
  ID = 1:5,
  Category = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")
)

# writeRaster(r_class,
#             filename = "08_Outputs/DyeCreek_S2_change_quintiles_pre_fire.tif",
#             filetype = "GTiff",  # GeoTIFF format
#             overwrite = TRUE)
# 
# writeRaster(r,
#             filename = "08_Outputs/DyeCreek_S2_change_raw_pre_fire.tif",
#             filetype = "GTiff",  # GeoTIFF format
#             overwrite = TRUE)

writeRaster(r,
            filename = "1_Data/08_Outputs/DyeCreek_S2_change_raw.tif",
            filetype = "GTiff",  # GeoTIFF format
            overwrite = TRUE)

writeRaster(r1,
            filename = "08_Outputs/DyeCreek_S2_change_raw_TMR1.tif",
            filetype = "GTiff",  # GeoTIFF format
            overwrite = TRUE)

