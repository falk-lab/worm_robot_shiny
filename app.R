# ================================
# SLC Lifespan Analysis (Merged + Plot Fixes)
# - Supports CSV or ZIP raw workflows
# - Robustly reads tab-delimited "csv" exports (like your reference file)
# - Uses cleaned_data as the shared source for all downstream analysis
# - Plots no longer silently disappear:
#     * validates data presence
#     * auto-uses available days if 2/5/10/15/20/25/30 aren't present
#     * shows clear messages instead of blank plot area
# ================================

library(shiny)
library(tidyverse)
library(vroom)
library(DT)
library(lubridate)
library(scales)
library(viridis)
library(ggh4x)
library(broom)
library(readxl)
library(platetools)
library(purrr)

ui <- fluidPage(
  titlePanel("SLC Lifespan Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("data_file", "Upload Data (.csv or tab-delimited .csv)", accept = c(".csv", ".tsv", ".txt")),
      fileInput("zip_file", "Upload ZIP Folder with Raw Data (.zip)", accept = ".zip"),
      fileInput("plate_map_file", "Upload Plate Map (.xlsx)", accept = ".xlsx"),
      fileInput("meta_file", "Upload Metadata (.csv)", accept = ".csv"),
      radioButtons(
        "data_source", "Data Source",
        choices = c(
          "Single CSV (Original)" = "csv",
          "ZIP Folder (Raw Robot Files)" = "zip"
        ),
        selected = "csv"
      ),
      textInput("ref_group", "Reference Group"),
      actionButton("run_btn", "Run Analysis", class = "btn-primary"),
      hr(),
      h4("Select Plot:"),
      selectInput(
        "plot_select", "Choose Plot",
        choices = c(
          "Raw Activity (Boxplot)" = "raw_box",
          "Activity (Worm Fraction Boxplot)" = "wf_box",
          "Activity (Worm Fraction Smoothed)" = "wf_smooth",
          "Z-score Activity (Boxplot)" = "z_box",
          "Z-score Activity (Smoothed)" = "z_smooth",
          "Day Max Activity" = "day_max",
          "Healthspan" = "healthspan",
          "Lifespan" = "lifespan",
          "Combined Plot (Day Max / Healthspan / Lifespan)" = "combined"
        )
      ),
      hr(),
      h4("Select Table:"),
      selectInput(
        "table_select", "Choose Table",
        choices = c(
          "Cleaned Data Table" = "cleaned_data",
          "Activity Ranked Table" = "activity_table",
          "Lifespan/Healthspan Table" = "lifespan_table",
          "Activity Worm Fraction Stats" = "activity_wf_stats",
          "Activity Z-score Stats" = "activity_zscore_stats",
          "Day Max Activity Table" = "day_max_table"
        )
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Plots", plotOutput("plot_main", height = "700px")),
        tabPanel("Tables", DTOutput("table_main"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  processed_data <- reactiveValues(
    cleaned_data = NULL,
    activity_ranked = NULL,
    lifespan = NULL,
    day_max_act = NULL,
    activity_wf_stats = NULL,
    activity_zscore_stats = NULL
  )
  
  observeEvent(input$run_btn, {
    req(input$data_source)
    
    # =============================
    # CSV workflow (robust to tab-delimited files)
    # =============================
    if (input$data_source == "csv") {
      req(input$data_file)
      
      # 1) Read attempt (comma)
      raw_data <- suppressMessages(vroom(input$data_file$datapath, show_col_types = FALSE))
      
      # 2) If it came in as a single column, re-read as tab-delimited
      if (ncol(raw_data) == 1) {
        raw_data <- suppressMessages(vroom(input$data_file$datapath, delim = "\t", show_col_types = FALSE))
      }
      
      # 3) Robust timestamp parsing to ensure days computes correctly
      raw_data <- raw_data %>%
        mutate(
          snapshot = as.character(snapshot),
          start_date = as.character(start_date),
          
          # handle ISO "2025-06-12T18:06:41Z" + variants
          snapshot = str_replace(snapshot, "T", " "),
          snapshot = str_remove(snapshot, "Z$"),
          snapshot = suppressWarnings(ymd_hms(snapshot, tz = "UTC")),
          snapshot = if_else(
            is.na(snapshot),
            suppressWarnings(parse_date_time(as.character(snapshot), orders = c("ymd HMS", "ymd HMSz"), tz = "UTC")),
            snapshot
          ),
          
          start_date = suppressWarnings(ymd(start_date, tz = "UTC")),
          
          picture = suppressWarnings(as.integer(picture)),
          well_num = suppressWarnings(as.integer(well_num)),
          worm_fraction = suppressWarnings(as.numeric(worm_fraction)),
          roi_activity = suppressWarnings(as.numeric(roi_activity)),
          plate = as.character(plate),
          well = as.character(well),
          well_info = as.character(well_info)
        )
      
      processed_data$cleaned_data <- raw_data
    }
    
    # =============================
    # ZIP workflow
    # =============================
    if (input$data_source == "zip") {
      req(input$zip_file)
      req(input$plate_map_file)
      
      # Well numbers mapping (A1-D6 -> 1-24)
      well24_to_num <- tibble(well = c(paste0('A', 1:6),
                                       paste0('B', 1:6),
                                       paste0('C', 1:6),
                                       paste0('D', 1:6))) %>%
        separate(well, into = c('well_row', 'well_col'), sep = 1, remove = FALSE) %>%
        arrange(well_col, well_row) %>%
        mutate(well_num = 1:24) %>%
        select(well, well_num)
      
      # Plate map parsing
      plate_maps <- readxl::read_excel(input$plate_map_file$datapath, sheet = 1, skip = 12) %>%
        rename(well_row = `...1`) %>%
        pivot_longer(`1`:`6`, names_to = 'well_col', values_to = 'well_info') %>%
        unite(well, c(well_row, well_col), sep = '') %>%
        select(well_info, well) %>%
        nest() %>%
        mutate(
          data = map(data, ~ platetools::plate_matrix(data = .$well_info, well = .$well, plate = 24)),
          data = map(data, ~ platetools::rotate_plate(.)),
          data = map(data, ~ as_tibble(.))
        ) %>%
        unnest(c(data)) %>%
        mutate(well_row = rep(LETTERS[1:4], nrow(.) / 4)) %>%
        pivot_longer(V1:V6, names_to = 'well_col', values_to = 'well_info') %>%
        mutate(well_col = str_remove(well_col, 'V')) %>%
        unite(well, c(well_row, well_col), sep = '') %>%
        left_join(well24_to_num, by = 'well') %>%
        mutate(plate = 1)
      
      # List robot files inside zip
      worm_robot_files <- unzip(input$zip_file$datapath, list = TRUE) %>%
        filter(str_detect(Name, 'Session')) %>%
        pull(Name)
      
      # Find measurement block locations
      measurement_locations <- NULL
      for (j in worm_robot_files) {
        temp <- vroom(unz(input$zip_file$datapath, j), delim = ',', show_col_types = FALSE) %>%
          mutate(start_rownum = row_number(),
                 file_path = j) %>%
          filter(str_detect(`<ROI Number>`, '<')) %>%
          mutate(measurement_clean = tolower(str_replace_all(str_remove_all(`<ROI Number>`, '<|>|,'), ' ', '_'))) %>%
          select(file_path, measurement_clean, measurement_name = `<ROI Number>`, start_rownum)
        
        temp <- temp %>%
          mutate(
            end_rownum = lead(start_rownum) - 1,
            date = str_remove_all(str_extract(file_path, '_[0-9]{4}-[0-9]{2}-[0-9]{2} '), '_| '),
            time = str_remove_all(str_extract(file_path, '\\([0-9]{2}-[0-9]{2}-[0-9]{2}\\)'), '\\(|\\)'),
            time = str_replace_all(time, '-', ':'),
            snapshot = paste0(date, ' ', time)
          )
        
        measurement_locations <- bind_rows(measurement_locations, temp)
      }
      
      # Worm fractions
      worm_frac_locations <- measurement_locations %>% filter(measurement_clean == 'worm_fraction')
      worm_frac <- NULL
      for (k in worm_robot_files) {
        temp_locations_worm_frac <- worm_frac_locations %>%
          filter(str_detect(file_path, fixed(k)))
        
        temp <- vroom(
          unz(input$zip_file$datapath, k),
          delim = ',',
          skip = temp_locations_worm_frac$start_rownum + 1,
          col_names = paste0(1:24),
          col_types = 'dddddddddddddddddddddddd',
          show_col_types = FALSE
        ) %>%
          head(1) %>%
          pivot_longer(
            `1`:`24`,
            names_to = 'well_num',
            values_to = 'worm_fraction',
            names_transform = list(well_num = as.integer)
          ) %>%
          mutate(snapshot = temp_locations_worm_frac$snapshot[1],
                 file_path = k)
        
        worm_frac <- bind_rows(worm_frac, temp)
      }
      
      # ROI activity
      activity_locations <- measurement_locations %>% filter(measurement_clean == 'roi_activity')
      roi_activity <- NULL
      for (l in worm_robot_files) {
        temp_locations_roi_activity <- activity_locations %>%
          filter(str_detect(file_path, fixed(l)))
        
        n_rows <- temp_locations_roi_activity$end_rownum[1] - temp_locations_roi_activity$start_rownum[1]
        
        temp <- vroom(
          unz(input$zip_file$datapath, l),
          delim = ',',
          skip = temp_locations_roi_activity$start_rownum[1] + 1,
          col_names = paste0(1:24),
          col_types = 'dddddddddddddddddddddddd',
          show_col_types = FALSE
        ) %>%
          head(n_rows) %>%
          mutate(picture = 1:n(),
                 file_path = l) %>%
          pivot_longer(
            `1`:`24`,
            names_to = 'well_num',
            values_to = 'roi_activity',
            names_transform = list(well_num = as.integer)
          ) %>%
          mutate(snapshot = temp_locations_roi_activity$snapshot[1])
        
        roi_activity <- bind_rows(roi_activity, temp)
      }
      
      # Combine into cleaned_data
      raw_data <- roi_activity %>%
        select(-file_path) %>%
        left_join(worm_frac %>% select(-file_path), by = join_by(snapshot, well_num)) %>%
        mutate(
          start_date = str_extract(snapshot, '[0-9]{4}-[0-9]{2}-[0-9]{2}'),
          start_date = ymd(start_date, tz = "UTC")
        ) %>%
        left_join(plate_maps, by = join_by(well_num)) %>%
        mutate(
          plate = if_else(is.na(plate), 1, plate),
          norm_activity = roi_activity * worm_fraction,
          snapshot = ymd_hms(snapshot, tz = "UTC")
        )
      
      processed_data$cleaned_data <- raw_data
    }
    
    # =============================
    # Downstream analysis (shared)
    # Uses processed_data$cleaned_data for BOTH workflows
    # =============================
    req(processed_data$cleaned_data)
    raw_data <- processed_data$cleaned_data
    
    activity_ranked <- raw_data %>%
      filter(!(picture %in% c(19, 31))) %>%
      group_by(plate, start_date, well, well_num, well_info, snapshot, worm_fraction) %>%
      summarize(mean_activity = mean(roi_activity, na.rm = TRUE), .groups = 'drop') %>%
      group_by(plate, start_date, well, well_num, well_info) %>%
      mutate(rescale_tot = scales::rescale(mean_activity)) %>%
      ungroup() %>%
      mutate(
        rescale_tot_inverse = 1 - rescale_tot,
        snapshot_date = as_date(snapshot),
        days = as.integer(snapshot_date - start_date)
      ) %>%
      filter(!is.na(well_info)) %>%
      separate(well_info, into = c('well_info2', 'replicate'),
               sep = "\\s+", remove = FALSE, fill = "right", extra = "merge") %>%
      mutate(
        replicate = str_remove_all(replace_na(replicate, ""), '\\(|\\)'),
        well_info2 = factor(well_info2, levels = c('N2', 'deletion', 'insertion', 'P107L', 'T239N')),
        activity_worm_frac = mean_activity * worm_fraction
      )
    
    lifespan <- activity_ranked %>%
      group_by(plate, start_date, well, well_num, well_info, well_info2, replicate) %>%
      summarise(
        t99_rank = nth(rescale_tot_inverse, which.min(abs(rescale_tot_inverse - 0.99))),
        t82_rank = nth(rescale_tot_inverse, which.min(abs(rescale_tot_inverse - 0.82))),
        .groups = 'drop'
      ) %>%
      left_join(
        select(activity_ranked, plate, start_date, well, well_num, well_info, days, rescale_tot_inverse),
        by = join_by(plate, start_date, well, well_num, well_info, t99_rank == rescale_tot_inverse)
      ) %>%
      rename(lifespan = days) %>%
      left_join(
        select(activity_ranked, plate, start_date, well, well_num, well_info, rescale_tot_inverse, days),
        by = join_by(plate, start_date, well, well_num, well_info, t82_rank == rescale_tot_inverse)
      ) %>%
      rename(healthspan = days) %>%
      group_by(plate, start_date, well, well_num, well_info) %>%
      arrange(lifespan) %>%
      filter(row_number() == 1) %>%
      ungroup()
    
    day_max_act <- activity_ranked %>%
      group_by(plate, start_date, well, well_num, well_info, well_info2, replicate, worm_fraction, days) %>%
      summarize(activity_worm_frac = mean_activity * worm_fraction, .groups = 'drop') %>%
      filter(!is.infinite(activity_worm_frac), !is.na(well)) %>%
      group_by(plate, well, well_info) %>%
      filter(activity_worm_frac == max(activity_worm_frac, na.rm = TRUE)) %>%
      ungroup()
    
    activity_wf_stats <- activity_ranked %>%
      group_by(plate, days) %>%
      nest() %>%
      mutate(test = map(data, ~ TukeyHSD(aov(activity_worm_frac ~ well_info, data = .)) %>% broom::tidy())) %>%
      select(-data) %>%
      unnest(test) %>%
      filter(adj.p.value < 0.05) %>%
      mutate(measurement = 'activity_wormfrac_norm') %>%
      select(measurement, plate, days, contrast, diff_means = estimate, adj.p.value)
    
    activity_zscore_stats <- activity_ranked %>%
      group_by(plate) %>%
      mutate(scaled_activity = as.numeric(scale(mean_activity))) %>%
      ungroup() %>%
      group_by(plate, days) %>%
      nest() %>%
      mutate(test = map(data, ~ TukeyHSD(aov(scaled_activity ~ well_info, data = .)) %>% broom::tidy())) %>%
      select(-data) %>%
      unnest(test) %>%
      filter(adj.p.value < 0.05) %>%
      mutate(measurement = 'activity_zscore') %>%
      select(measurement, plate, days, contrast, diff_means = estimate, adj.p.value)
    
    processed_data$activity_ranked <- activity_ranked
    processed_data$lifespan <- lifespan
    processed_data$day_max_act <- day_max_act
    processed_data$activity_wf_stats <- activity_wf_stats
    processed_data$activity_zscore_stats <- activity_zscore_stats
  })
  
  # --- Plot Rendering (fixed) ---
  output$plot_main <- renderPlot({
    req(processed_data$activity_ranked)
    
    activity_ranked <- processed_data$activity_ranked
    lifespan <- processed_data$lifespan
    day_max_act <- processed_data$day_max_act
    plot_sel <- input$plot_select
    
    validate(
      need(nrow(activity_ranked) > 0, "No activity data available after processing."),
      need(!all(is.na(activity_ranked$days)), "Days could not be computed (snapshot/start_date parsing issue).")
    )
    
    preferred_days <- c(2, 5, 10, 15, 20, 25, 30)
    available_days <- sort(unique(activity_ranked$days[!is.na(activity_ranked$days)]))
    days_to_use <- intersect(preferred_days, available_days)
    if (length(days_to_use) == 0) days_to_use <- head(available_days, 7)
    
    switch(plot_sel,
           "raw_box" = {
             df <- activity_ranked %>% filter(days %in% days_to_use)
             validate(need(nrow(df) > 0, "No rows available for Raw Activity plot."))
             ggplot(df, aes(x = well_info2, y = mean_activity, color = well_info2)) +
               geom_boxplot() +
               facet_grid2(plate ~ days, scales = "free_x", independent = "x") +
               theme_bw() +
               labs(x = "Condition", y = "Raw Activity", color = "Condition")
           },
           "wf_box" = {
             df <- activity_ranked %>% filter(days %in% days_to_use) %>%
               mutate(activity_worm_frac = mean_activity * worm_fraction)
             validate(need(nrow(df) > 0, "No rows available for Worm Fraction boxplot."))
             ggplot(df, aes(x = well_info2, y = activity_worm_frac, color = well_info2)) +
               geom_boxplot() +
               facet_grid2(plate ~ days, scales = "free_x", independent = "x") +
               theme_bw() +
               labs(x = "Condition", y = "Activity × Worm Fraction")
           },
           "wf_smooth" = {
             df <- activity_ranked %>% mutate(activity_worm_frac = mean_activity * worm_fraction)
             validate(need(nrow(df) > 0, "No rows available for Worm Fraction smoothed plot."))
             ggplot(df, aes(x = days, y = activity_worm_frac, color = well_info2, fill = well_info2)) +
               geom_smooth(alpha = 0.2) +
               facet_wrap(~ plate) +
               theme_bw() +
               labs(x = "Days", y = "Activity × Worm Fraction", color = "Condition", fill = "Condition")
           },
           "z_box" = {
             df <- activity_ranked %>%
               group_by(plate) %>%
               mutate(scaled_activity = as.numeric(scale(mean_activity))) %>%
               ungroup() %>%
               filter(days %in% days_to_use)
             validate(need(nrow(df) > 0, "No rows available for Z-score boxplot."))
             ggplot(df, aes(x = well_info2, y = scaled_activity, color = well_info2)) +
               geom_boxplot() +
               facet_grid2(plate ~ days, scales = "free_x", independent = "x") +
               theme_bw() +
               labs(x = "Condition", y = "Z-score Activity", color = "Condition")
           },
           "z_smooth" = {
             df <- activity_ranked %>%
               group_by(plate) %>%
               mutate(scaled_activity = as.numeric(scale(mean_activity))) %>%
               ungroup()
             validate(need(nrow(df) > 0, "No rows available for Z-score smoothed plot."))
             ggplot(df, aes(x = days, y = scaled_activity, color = well_info2, fill = well_info2)) +
               geom_smooth(alpha = 0.2) +
               facet_wrap(~ plate) +
               theme_bw() +
               labs(x = "Days", y = "Z-score Activity", color = "Condition", fill = "Condition")
           },
           "day_max" = {
             validate(need(!is.null(day_max_act) && nrow(day_max_act) > 0, "Day max activity table is empty."))
             ggplot(day_max_act, aes(x = well_info2, y = days)) +
               geom_boxplot() +
               theme_bw() +
               labs(x = "Condition", y = "Day of Max Activity")
           },
           "healthspan" = {
             validate(need(!is.null(lifespan) && nrow(lifespan) > 0, "Lifespan table is empty."))
             ggplot(lifespan, aes(x = well_info2, y = healthspan)) +
               geom_boxplot() +
               theme_bw() +
               labs(x = "Condition", y = "Healthspan (Days)")
           },
           "lifespan" = {
             validate(need(!is.null(lifespan) && nrow(lifespan) > 0, "Lifespan table is empty."))
             ggplot(lifespan, aes(x = well_info2, y = lifespan)) +
               geom_boxplot() +
               theme_bw() +
               labs(x = "Condition", y = "Lifespan (Days)")
           },
           "combined" = {
             validate(need(!is.null(lifespan) && nrow(lifespan) > 0, "Lifespan table is empty."))
             df <- lifespan %>%
               left_join(
                 day_max_act %>% select(plate, well, well_info, well_info2, days_max = days),
                 by = c("plate", "well", "well_info", "well_info2")
               ) %>%
               pivot_longer(cols = c(days_max, healthspan, lifespan),
                            names_to = "measurement", values_to = "days") %>%
               mutate(measurement = recode(measurement,
                                           days_max = "Day Max Activity",
                                           healthspan = "Healthspan",
                                           lifespan = "Lifespan"))
             validate(need(nrow(df) > 0, "No rows available for combined plot."))
             ggplot(df, aes(x = well_info2, y = days, color = measurement)) +
               geom_boxplot() +
               scale_color_viridis_d(option = "rocket", end = 0.75) +
               theme_bw() +
               labs(x = "Condition", y = "Days", color = "")
           }
    )
  })
  
  # --- Table Rendering ---
  output$table_main <- renderDT({
    req(processed_data$cleaned_data)
    
    tbl_sel <- input$table_select
    ref_group <- trimws(input$ref_group)
    
    switch(tbl_sel,
           "cleaned_data" = {
             df <- processed_data$cleaned_data
             if (!is.null(ref_group) && ref_group != "") df <- df %>% filter(well_info == ref_group)
             datatable(df, extensions = 'Buttons', options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel')))
           },
           "activity_table" = {
             req(processed_data$activity_ranked)
             df <- processed_data$activity_ranked %>%
               select(plate, start_date, well, well_info, days, snapshot, worm_fraction, mean_activity, activity_worm_frac)
             if (!is.null(ref_group) && ref_group != "") df <- df %>% filter(well_info == ref_group)
             datatable(df, extensions = 'Buttons', options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel')))
           },
           "lifespan_table" = {
             req(processed_data$lifespan)
             df <- processed_data$lifespan %>% select(plate, start_date, well, well_info, healthspan, lifespan)
             if (!is.null(ref_group) && ref_group != "") df <- df %>% filter(well_info == ref_group)
             datatable(df, extensions = 'Buttons', options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel')))
           },
           "day_max_table" = {
             req(processed_data$day_max_act)
             df <- processed_data$day_max_act %>% select(plate, start_date, well, well_info, days)
             if (!is.null(ref_group) && ref_group != "") df <- df %>% filter(well_info == ref_group)
             datatable(df, extensions = 'Buttons', options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel')))
           },
           "activity_wf_stats" = {
             req(processed_data$activity_wf_stats)
             datatable(processed_data$activity_wf_stats, extensions = 'Buttons',
                       options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel')))
           },
           "activity_zscore_stats" = {
             req(processed_data$activity_zscore_stats)
             datatable(processed_data$activity_zscore_stats, extensions = 'Buttons',
                       options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel')))
           }
    )
  })
}

shinyApp(ui = ui, server = server)