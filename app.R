library(shiny)
library(tidyverse)
library(vroom)
library(DT)
library(lubridate)
library(scales)
library(viridis)
library(ggh4x)
library(broom)

ui <- fluidPage(
  titlePanel("SLC Lifespan Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("data_file", "Upload Data (.csv)", accept = ".csv"),
      fileInput("meta_file", "Upload Metadata (.csv)", accept = ".csv"),
      textInput("ref_group", "Reference Group"),
      actionButton("run_btn", "Run Analysis", class = "btn-primary"),
      hr(),
      h4("Select Plot:"),
      selectInput("plot_select", "Choose Plot",
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
                  )),
      hr(),
      h4("Select Table:"),
      selectInput("table_select", "Choose Table",
                  choices = c(
                    "Activity Ranked Table" = "activity_table",
                    "Lifespan/Healthspan Table" = "lifespan_table",
                    "Activity Worm Fraction Stats" = "activity_wf_stats",
                    "Activity Z-score Stats" = "activity_zscore_stats",
                    "Day Max Activity Table" = "day_max_table"
                  ))
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
    activity_ranked = NULL,
    lifespan = NULL,
    day_max_act = NULL,
    activity_wf_stats = NULL,
    activity_zscore_stats = NULL
  )
  
  observeEvent(input$run_btn, {
    req(input$data_file)
    
    raw_data <- vroom(input$data_file$datapath)
    
    # --- Activity Processing ---
    activity_ranked <- raw_data %>%
      filter(!(picture %in% c(19, 31))) %>% 
      group_by(plate, start_date, well, well_num, well_info, snapshot, worm_fraction) %>%
      summarize(mean_activity = mean(roi_activity), .groups = 'drop') %>%
      group_by(plate, start_date, well, well_num, well_info) %>%
      mutate(rescale_tot = scales::rescale(mean_activity)) %>%
      ungroup() %>%
      mutate(rescale_tot_inverse = 1 - rescale_tot,
             snapshot_date = lubridate::as_date(snapshot),
             days = as.integer(snapshot_date - start_date)) %>%
      filter(!is.na(well_info)) %>%
      separate(well_info, into = c('well_info2', 'replicate'), sep = ' ', remove = FALSE) %>%
      mutate(replicate = str_remove_all(replicate, '\\(|\\)'),
             well_info2 = factor(well_info2, levels = c('N2', 'deletion', 'insertion',
                                                        'P107L', 'T239N')),
             activity_worm_frac = mean_activity * worm_fraction)
    
    # --- Lifespan/Healthspan ---
    lifespan <- activity_ranked %>%
      group_by(plate, start_date, well, well_num, well_info, well_info2, replicate) %>%
      summarise(t99_rank = nth(rescale_tot_inverse, which.min(abs(rescale_tot_inverse-0.99))),
                t82_rank = nth(rescale_tot_inverse, which.min(abs(rescale_tot_inverse-0.82))),
                .groups = 'drop') %>%
      left_join(select(activity_ranked, plate, start_date, well, well_num, well_info,
                       days, rescale_tot_inverse),
                by = join_by(plate, start_date, well, well_num, well_info,
                             t99_rank == rescale_tot_inverse)) %>%
      rename(lifespan = days) %>%
      left_join(select(activity_ranked, plate, start_date, well, well_num, well_info,
                       rescale_tot_inverse, days),
                by = join_by(plate, start_date, well, well_num, well_info,
                             t82_rank == rescale_tot_inverse)) %>%
      rename(healthspan = days) %>%
      group_by(plate, start_date, well, well_num, well_info) %>%
      arrange(lifespan) %>%
      filter(row_number() == 1) %>%
      ungroup()
    
    # --- Day Max Activity ---
    day_max_act <- activity_ranked %>%
      group_by(plate, start_date, well, well_num, well_info, well_info2, replicate,
               worm_fraction, days) %>%
      summarize(activity_worm_frac = mean_activity * worm_fraction) %>%
      ungroup() %>%
      filter(!is.infinite(activity_worm_frac), !is.na(well)) %>%
      group_by(plate, well, well_info) %>%
      filter(activity_worm_frac == max(activity_worm_frac)) %>%
      ungroup()
    
    # --- WF Stats (Tukey HSD) ---
    activity_wf_stats <- activity_ranked %>%
      group_by(plate, days) %>%
      nest() %>%
      mutate(test = map(data, ~ TukeyHSD(aov(activity_worm_frac ~ well_info, data = .)) %>%
                          broom::tidy())) %>%
      select(-data) %>%
      unnest(test) %>%
      filter(adj.p.value < 0.05) %>%
      mutate(measurement = 'activity_wormfrac_norm') %>%
      select(measurement, plate, days, contrast, diff_means = estimate, adj.p.value)
    
    # --- Z-score Stats (Tukey HSD) ---
    activity_zscore_stats <- activity_ranked %>%
      group_by(plate) %>%
      mutate(scaled_activity = scale(mean_activity)) %>%
      ungroup() %>%
      group_by(plate, days) %>%
      nest() %>%
      mutate(test = map(data, ~ TukeyHSD(aov(scaled_activity ~ well_info, data = .)) %>%
                          broom::tidy())) %>%
      select(-data) %>%
      unnest(test) %>%
      filter(adj.p.value < 0.05) %>%
      mutate(measurement = 'activity_zscore') %>%
      select(measurement, plate, days, contrast, diff_means = estimate, adj.p.value)
    
    # Store data
    processed_data$activity_ranked <- activity_ranked
    processed_data$lifespan <- lifespan
    processed_data$day_max_act <- day_max_act
    processed_data$activity_wf_stats <- activity_wf_stats
    processed_data$activity_zscore_stats <- activity_zscore_stats
  })
  
  # --- Plot Rendering ---
  output$plot_main <- renderPlot({
    req(processed_data$activity_ranked)
    
    activity_ranked <- processed_data$activity_ranked
    lifespan <- processed_data$lifespan
    day_max_act <- processed_data$day_max_act
    plot_sel <- input$plot_select
    
    switch(plot_sel,
           "raw_box" = {
             ggplot(activity_ranked %>% filter(days %in% c(2,5,10,15,20,25,30)),
                    aes(x = well_info2, y = mean_activity, color = well_info2)) +
               geom_boxplot() +
               facet_grid2(plate ~ days, scales = "free_x", independent = 'x') +
               theme_bw() +
               labs(x = 'Condition', y = 'Raw Activity', color = 'Condition')
           },
           "wf_box" = {
             ggplot(activity_ranked %>%
                      filter(days %in% c(2,5,10,15,20,25,30)) %>%
                      mutate(activity_worm_frac = mean_activity * worm_fraction),
                    aes(x = well_info2, y = activity_worm_frac, color = well_info2)) +
               geom_boxplot() +
               facet_grid2(plate ~ days, scales = "free_x", independent = 'x') +
               theme_bw() +
               labs(x = 'Condition', y = 'Activity × Worm Fraction')
           },
           "wf_smooth" = {
             ggplot(activity_ranked %>%
                      mutate(activity_worm_frac = mean_activity * worm_fraction),
                    aes(x = days, y = activity_worm_frac, color = well_info2, fill = well_info2)) +
               geom_smooth(alpha = 0.2) +
               facet_wrap(~ plate) +
               theme_bw() +
               labs(x = 'Days', y = 'Activity × Worm Fraction', color = 'Condition', fill = 'Condition')
           },
           "z_box" = {
             ggplot(activity_ranked %>%
                      group_by(plate) %>%
                      mutate(scaled_activity = scale(mean_activity)) %>%
                      ungroup() %>%
                      filter(days %in% c(2,5,10,15,20,25,30)),
                    aes(x = well_info2, y = scaled_activity, color = well_info2)) +
               geom_boxplot() +
               facet_grid2(plate ~ days, scales = "free_x", independent = 'x') +
               theme_bw() +
               labs(x = 'Condition', y = 'Z-score Activity', color = 'Condition')
           },
           "z_smooth" = {
             ggplot(activity_ranked %>%
                      group_by(plate) %>%
                      mutate(scaled_activity = scale(mean_activity)) %>%
                      ungroup(),
                    aes(x = days, y = scaled_activity, color = well_info2, fill = well_info2)) +
               geom_smooth(alpha = 0.2) +
               facet_wrap(~ plate) +
               theme_bw() +
               labs(x = 'Days', y = 'Z-score Activity', color = 'Condition', fill = 'Condition')
           },
           "day_max" = {
             ggplot(day_max_act, aes(x = well_info2, y = days)) +
               geom_boxplot() +
               theme_bw() +
               labs(x = 'Condition', y = 'Day of Max Activity')
           },
           "healthspan" = {
             ggplot(lifespan, aes(x = well_info2, y = healthspan)) +
               geom_boxplot() +
               theme_bw() +
               labs(x = 'Condition', y = 'Healthspan (Days)')
           },
           "lifespan" = {
             ggplot(lifespan, aes(x = well_info2, y = lifespan)) +
               geom_boxplot() +
               theme_bw() +
               labs(x = 'Condition', y = 'Lifespan (Days)')
           },
           "combined" = {
             lifespan %>%
               left_join(day_max_act %>%
                           select(plate, well, well_info, well_info2, days_max = days),
                         by = c("plate","well","well_info","well_info2")) %>%
               pivot_longer(cols = c(days_max, healthspan, lifespan),
                            names_to = "measurement", values_to = "days") %>%
               mutate(measurement = recode(measurement,
                                           days_max = "Day Max Activity",
                                           healthspan = "Healthspan",
                                           lifespan = "Lifespan")) %>%
               ggplot(aes(x = well_info2, y = days, color = measurement)) +
               geom_boxplot() +
               scale_color_viridis_d(option = "rocket", end = 0.75) +
               theme_bw() +
               labs(x = 'Condition', y = 'Days', color = '')
           })
  })
  
  # --- Tables (unchanged) ---
  output$table_main <- renderDT({
    req(processed_data$activity_ranked)
    table_sel <- input$table_select
    ref_group <- trimws(input$ref_group)
    
    switch(table_sel,
           "activity_table" = {
             df <- processed_data$activity_ranked %>%
               select(plate:well, well_info, days, snapshot, worm_fraction,
                      mean_activity, activity_worm_frac)
             if (ref_group != "") df <- df %>% filter(well_info == ref_group)
             datatable(df, extensions='Buttons',
                       options=list(dom='Blfrtip', buttons=c('copy','csv','excel')))
           },
           "lifespan_table" = {
             df <- processed_data$lifespan
             if (ref_group != "") df <- df %>% filter(well_info == ref_group)
             datatable(df, extensions='Buttons',
                       options=list(dom='Blfrtip', buttons=c('copy','csv','excel')))
           },
           "activity_wf_stats" = datatable(processed_data$activity_wf_stats,
                                           extensions='Buttons',
                                           options=list(dom='Blfrtip', buttons=c('copy','csv','excel'))),
           "activity_zscore_stats" = datatable(processed_data$activity_zscore_stats,
                                               extensions='Buttons',
                                               options=list(dom='Blfrtip', buttons=c('copy','csv','excel'))),
           "day_max_table" = datatable(processed_data$day_max_act,
                                       extensions='Buttons',
                                       options=list(dom='Blfrtip', buttons=c('copy','csv','excel')))
    )
  })
}

shinyApp(ui = ui, server = server)
