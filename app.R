library(shiny)
library(tidyverse)
library(vroom)
library(DT)
library(lubridate)
library(scales)
library(viridis)
library(ggh4x)

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
                    "Lifespan/Healthspan Table" = "lifespan_table"
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
  
  # Reactive value to store processed data
  processed_data <- reactiveValues(activity_ranked = NULL, lifespan = NULL, day_max_act = NULL)
  
  observeEvent(input$run_btn, {
    req(input$data_file)
    
    # Read in data
    raw_data <- vroom(input$data_file$datapath)
    
    # Activity Data Processing
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
    
    # Lifespan / Healthspan
    lifespan <- activity_ranked %>%
      group_by(plate, start_date, well, well_num, well_info, well_info2, replicate) %>%
      summarise(t99_rank = nth(rescale_tot_inverse, 
                               which.min(abs(rescale_tot_inverse-0.99))),
                t82_rank = nth(rescale_tot_inverse, 
                               which.min(abs(rescale_tot_inverse-0.82))), .groups = 'drop') %>%
      left_join(select(activity_ranked,
                       plate, start_date, well, well_num, well_info, days,
                       rescale_tot_inverse),
                by = join_by(plate, start_date, well, well_num, well_info,
                             t99_rank == rescale_tot_inverse)) %>%
      rename(lifespan = days) %>%
      left_join(select(activity_ranked,
                       plate, start_date, well, well_num, well_info,
                       rescale_tot_inverse, days),
                by = join_by(plate, start_date, well, well_num, well_info,
                             t82_rank == rescale_tot_inverse)) %>%
      rename(healthspan = days) %>%
      group_by(plate, start_date, well, well_num, well_info) %>%
      arrange(lifespan) %>%
      filter(row_number() == 1) %>%
      ungroup()
    
    # Day Max Activity
    day_max_act <- activity_ranked %>%
      group_by(plate, start_date, well, well_num, well_info, well_info2, replicate,
               worm_fraction, days) %>%
      summarize(activity_worm_frac = mean_activity * worm_fraction, .groups = 'drop') %>%
      filter(!is.infinite(activity_worm_frac), !is.na(well)) %>%
      group_by(plate, well, well_info) %>%
      filter(activity_worm_frac == max(activity_worm_frac)) %>%
      ungroup()
    
    # Store processed data
    processed_data$activity_ranked <- activity_ranked
    processed_data$lifespan <- lifespan
    processed_data$day_max_act <- day_max_act
  })
  
  # Render plot based on selection
  output$plot_main <- renderPlot({
    req(processed_data$activity_ranked)
    
    activity_ranked <- processed_data$activity_ranked
    lifespan <- processed_data$lifespan
    day_max_act <- processed_data$day_max_act
    
    plot_sel <- input$plot_select
    
    switch(plot_sel,
           "raw_box" = {
             activity_ranked %>% 
               filter(days %in% c(2,5,10,15,20,25,30)) %>%
               ggplot(aes(x = well_info2, y = mean_activity, color = well_info2)) +
               geom_boxplot() +
               ggh4x::facet_grid2(plate ~ days, scales = "free_x", independent = 'x') +
               labs(x = 'Experimental Condition',
                    y = 'Raw Activity',
                    color = 'condition',
                    title = 'Raw Activity') +
               theme_bw()
           },
           "wf_box" = {
             activity_ranked %>% 
               group_by(plate, start_date, well, well_num, well_info, well_info2, worm_fraction, days) %>%
               summarize(activity_worm_frac = mean_activity * worm_fraction, .groups = 'drop') %>%
               filter(days %in% c(2,5,10,15,20,25,30)) %>%
               ggplot(aes(x = well_info2, y = activity_worm_frac, color = well_info2)) +
               geom_boxplot() +
               ggh4x::facet_grid2(plate ~ days, scales = "free_x", independent = 'x') +
               labs(x = 'RNAi Concentration, Amino Acid Treatment',
                    y = 'Activity (Normalized to Worm Fraction)',
                    color = 'condition',
                    title = 'Worm Fraction Normalized Activity') +
               theme_bw()
           },
           "wf_smooth" = {
             activity_ranked %>% 
               group_by(plate, start_date, well, well_num, well_info, well_info2, worm_fraction, days) %>%
               summarize(activity_worm_frac = mean_activity * worm_fraction, .groups = 'drop') %>%
               ggplot(aes(x = days, y = activity_worm_frac, color = well_info2, fill = well_info2)) +
               geom_smooth(alpha = 0.2) +
               facet_wrap(~ plate) +
               labs(x = 'Time (Days)', y = 'Activity',
                    color = 'condition', fill = 'condition',
                    title = 'Worm Fraction Normalized Activity') +
               theme_bw()
           },
           "z_box" = {
             activity_ranked %>%
               group_by(plate) %>%
               mutate(scaled_activity = scale(mean_activity, center = TRUE, scale = TRUE)) %>%
               ungroup() %>%
               filter(days %in% c(2,5,10,15,20,25,30)) %>%
               ggplot(aes(x = well_info2, y = scaled_activity, color = well_info2)) +
               geom_boxplot() +
               ggh4x::facet_grid2(plate ~ days, scales = "free_x", independent = 'x') +
               labs(x = 'RNAi Concentration, Amino Acid Treatment',
                    y = 'Z-score Activity',
                    color = 'condition',
                    title = 'Z-score Normalized Activity') +
               theme_bw()
           },
           "z_smooth" = {
             activity_ranked %>%
               group_by(plate) %>%
               mutate(scaled_activity = scale(mean_activity, center = TRUE, scale = TRUE)) %>%
               ungroup() %>%
               ggplot(aes(x = days, y = scaled_activity, color = well_info2, fill = well_info2)) +
               geom_smooth(alpha = 0.2) +
               facet_wrap(~ plate) +
               labs(x = 'Time (Days)', y = 'Z-score Activity',
                    color = 'condition', fill = 'condition',
                    title = 'Z-score Normalized Activity') +
               theme_bw()
           },
           "day_max" = {
             day_max_act %>%
               ggplot(aes(x = well_info2, y = days)) +
               geom_boxplot() +
               coord_cartesian(ylim = c(0,10)) +
               labs(x = 'Mutant', y = 'Day of Max Activity') +
               theme_bw()
           },
           "healthspan" = {
             lifespan %>%
               ggplot(aes(x = well_info2, y = healthspan)) +
               geom_boxplot() +
               labs(x = 'Mutant', y = 'Healthspan (Days)', title = 'Healthspan') +
               theme_bw()
           },
           "lifespan" = {
             lifespan %>%
               ggplot(aes(x = well_info2, y = lifespan)) +
               geom_boxplot() +
               labs(x = 'Mutant', y = 'Lifespan (Days)', title = 'Lifespan') +
               theme_bw()
           },
           "combined" = {
             lifespan %>%
               left_join(day_max_act %>% select(plate, well, well_info, well_info2, days_max = days),
                         by = c("plate","well","well_info","well_info2")) %>%
               pivot_longer(cols = c(days_max, healthspan, lifespan),
                            names_to = "measurement", values_to = "days") %>%
               mutate(measurement = ifelse(measurement == "days_max","Day Max Activity",str_to_title(measurement))) %>%
               ggplot(aes(x = well_info2, y = days, color = measurement)) +
               geom_boxplot() +
               viridis::scale_color_viridis(discrete = TRUE, option = 'rocket', end = 0.75) +
               coord_cartesian(ylim = c(0,30)) +
               labs(x = 'RNAi Concentration, Amino Acid Treatment', y = 'Days', color = '') +
               theme_bw()
           }
    )
  })
  
  # Render table based on selection
  output$table_main <- renderDT({
    req(processed_data$activity_ranked)
    
    table_sel <- input$table_select
    
    switch(table_sel,
           "activity_table" = {
             processed_data$activity_ranked %>%
               select(plate:well, well_info, days, snapshot, worm_fraction, mean_activity,
                      activity_worm_frac) %>%
               datatable(extensions = 'Buttons', options = list(dom='Blfrtip',
                                                                buttons=c('copy','csv','excel'),
                                                                lengthMenu=list(c(10,25,50,-1),
                                                                                c(10,25,50,"All"))))
           },
           "lifespan_table" = {
             processed_data$lifespan %>%
               select(plate, start_date, well, well_info, healthspan, lifespan,
                      t82_rank, t99_rank) %>%
               datatable(extensions = 'Buttons', options = list(dom='Blfrtip',
                                                                buttons=c('copy','csv','excel'),
                                                                lengthMenu=list(c(10,25,50,-1),
                                                                                c(10,25,50,"All"))))
           }
    )
  })
}

shinyApp(ui = ui, server = server)
