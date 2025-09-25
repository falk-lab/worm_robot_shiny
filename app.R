library(shiny)

ui <- fluidPage()

server <- function(input, output, session) {
  
  # -------- Modal 1: Uploader --------
  uploader_modal <- modalDialog(
    easyClose = TRUE, footer = NULL,
    tags$div(
      style = "padding:12px 12px 24px; font-family:sans-serif; width:900px; max-width:95%; height:600px;",
      tags$h1("Start with Data Upload",
              style = "font-weight:700; font-size:48px; margin:0 0 16px 0;"),
      tags$div(
        style = "border:6px solid #000; border-radius:4px; height:520px; padding:24px;",
        tags$div(
          style = "border:4px solid #000; width:260px; padding:16px; border-radius:4px;",
          tags$h3("Upload Data:", style="margin-top:0;"),
          fileInput("data_file", label = NULL, accept = c(".csv",".tsv",".xlsx")),
          tags$h3("Upload Metadata:"),
          fileInput("meta_file", label = NULL, accept = c(".csv",".tsv",".xlsx")),
          downloadButton("dl_template", "Download template"),
          tags$h3("Reference group:", style="margin-top:18px;"),
          textInput("ref_group", label = NULL, placeholder = ""),
          actionButton("run_btn", "Run", class = "btn-primary")
        )
      )
    )
  )
  
  # -------- Modal 2: Analysis View --------
  analysis_modal <- modalDialog(
    easyClose = TRUE, footer = NULL,
    tags$div(
      style = "padding:12px 12px 24px; font-family:sans-serif; width:1100px; max-width:98%; height:700px;",
      tags$h1("Analysis View",
              style="font-weight:700; font-size:48px; margin:0 0 16px 0;"),
      tags$div(
        style="border:6px solid #000; border-radius:4px; height:580px; padding:24px; display:flex; gap:18px;",
        # Left panel
        tags$div(
          style="border:4px solid #000; width:300px; padding:18px; border-radius:4px; display:flex; flex-direction:column; gap:12px;",
          tags$h2("Upload Data:", style="margin:0; font-size:28px;"),
          fileInput("data_file2", label = NULL),
          tags$h2("Upload Metadata:", style="margin:0; font-size:28px;"),
          fileInput("meta_file2", label = NULL),
          downloadButton("dl_template2", "Download template"),
          tags$h2("Reference group:", style="margin:0; font-size:28px;"),
          textInput("ref_group2", label = NULL, placeholder = ""),
          actionButton("run_btn2", "Run"),
          tags$hr(),
          downloadButton("dl_report", "Download report")
        ),
        # Right area
        tags$div(
          style="flex:1; display:flex; flex-direction:column; gap:10px;",
          tags$div(
            style="display:flex; gap:12px;",
            actionButton("plot1_btn", "Plot 1"),
            actionButton("plot2_btn", "Plot 2"),
            actionButton("plot3_btn", "Plot 3"),
            actionButton("plot4_btn", "Plot 4")
          ),
          tags$div(
            style="border:4px solid #000; flex:1; display:flex; align-items:center; justify-content:center;",
            plotOutput("plot_main", height = "100%")
          )
        )
      )
    )
  )
  
  # Show uploader on load
  observe({ showModal(uploader_modal) })
  
  # Switch to analysis on Run
  observeEvent(input$run_btn, {
    removeModal()
    showModal(analysis_modal)
  })
  
  # Downloads
  output$dl_template <- downloadHandler(
    filename = function() "metadata_template.csv",
    content = function(file) {
      write.csv(data.frame(SampleID = character(0), Group = character(0)),
                file, row.names = FALSE)
    }
  )
  
  output$dl_template2 <- downloadHandler(
    filename = function() "metadata_template.csv",
    content = function(file) {
      write.csv(data.frame(SampleID = character(0), Group = character(0)),
                file, row.names = FALSE)
    }
  )
  
  output$dl_report <- downloadHandler(
    filename = function() "report.txt",
    content = function(file) {
      writeLines(c(
        "Example analysis report",
        paste("Reference group:",
              if (nzchar(input$ref_group)) input$ref_group else input$ref_group2),
        "…add your real results here…"
      ), file)
    }
  )
  
  # --- Plot switching logic ---
  current_plot <- reactiveVal("plot1")
  
  observeEvent(input$plot1_btn, current_plot("plot1"))
  observeEvent(input$plot2_btn, current_plot("plot2"))
  observeEvent(input$plot3_btn, current_plot("plot3"))
  observeEvent(input$plot4_btn, current_plot("plot4"))
  
  output$plot_main <- renderPlot({
    switch(current_plot(),
           plot1 = plot(cars, main = "Plot 1: cars dataset"),
           plot2 = plot(pressure, main = "Plot 2: pressure dataset"),
           plot3 = hist(cars$speed, main = "Plot 3: speed histogram"),
           plot4 = hist(pressure$temperature, main = "Plot 4: temp histogram")
    )
  })
}

shinyApp(ui = ui, server = server)
