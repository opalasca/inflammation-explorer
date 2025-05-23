
library("shiny")
library("DT")
library("readr")
library("shinyjs")
library("shinycssloaders")


options(shiny.maxRequestSize = 100 * 1024^2)

source("functions.R")


# Define UI ---------
ui <- fluidPage(
  
  useShinyjs(),  # Enable shinyjs in UI
  titlePanel("Explore inflammation in your dataset"),
  
  # Add info instructions right below the title
  tags$div(style = "margin-bottom: 20px; font-size: 14px; padding: 10px; background-color: #f9f9f9; border: 1px solid #ddd; border-radius: 5px;",
           tags$p("This app makes use of the gene sets derived in the paper 
                  “Mining inflammation in human transcriptomes: A consensus signature accross diseases and tissues” to check whether inflammation is 
                  impacting differential expression results."),
           
           # More Info button + hidden section (now directly after intro)
           actionButton("toggle_info", "More Info", icon = icon("info-circle"), class = "btn btn-outline-secondary btn-sm", style = "margin-bottom: 10px;"),
           
           shinyjs::hidden(
             div(id = "more_info_text", style = "margin-top: 10px;",
                 tags$p("We used rank aggregation to combine differential expression analysis results from 33 case–control comparisons spanning 18 inflammatory diseases and 9 different tissues. 
                        From this, we derived a high-confidence inflammation signature comprising the top 100 most consistently upregulated genes across all datasets. Whether this set shows positive enrichment can be tested with GSEA to check the presence of an inflammatory signal. We include a selection of 12 other inflammation-related terms from multiple ontology resources for comparison. 
                        We also defined a broader inflammatome, comprising the top 2000 most consistently upregulated genes. 
                        This set can be visualized in a volcano plot to assess the extent to which inflammation is driving gene or protein abundance differences in your dataset.")
             )
           ),
           tags$p(strong("How to use:")),
           tags$p("Upload your .txt, .tsv, .csv or .xlsx file with results from differential expression analysis. The file should contain at least the following columns: gene/protein identifier, log2 fold change, p-value or adjusted p-value (FDR)."),
           tags$p("Check the 'Data Preview' tab to help in selecting columns to be used in the analysis."),
           tags$p("• ", strong("ID column:"), " the column you want to use as gene identifier. It should be one of Ensembl Gene, Entrez, Uniprot, RefSeq or gene name (symbol)"),
           tags$p("• ", strong("Keytype:"), " the identifier type for the selected ID column. For example, if your file includes both UniprotID and gene symbol, and you select gene symbol as ID column, then keytype should be “Symbol”."),
           tags$p("• ", strong("Sorting column:"), " the column used for ranking genes in GSEA; preferably a test statistic like “stat” (DESeq2) or “t” (Limma), or logFC."),
           tags$p("• ", strong("logFC"), "(log fold change) and ", strong("P-value"), "(or adjusted p-value) columns for plotting.")
  ),
  
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload File", accept = c(".tsv", ".txt", ".csv",".xlsx")),
      # Dynamically show the sheet selection input if an Excel file is uploaded
      uiOutput("sheet_select_ui"),  # This is the key part to add the sheet select input dynamically
      
      div(style = "margin-top: 10px; font-weight: bold;", "...or use an example dataset:"),
      selectInput("example_file_choice", NULL, choices = 
                    c("Select Example Dataset" = "", 
                      "Park et. al, 2022, IgA nephropathy, RNA-Seq", 
                      "Bennike et. al, 2015, ulcerative colitis, proteomics")),
      
      # Add spacing before the next section
      hr(),
      
      selectInput("id_col", "Select ID Column", choices = NULL),
      selectInput("keytype", "Select Keytype", choices = c("Ensembl", "Uniprot", "Entrez", "Refseq", "Symbol")),
      selectInput("sort_col", "Select Sorting Column", choices = NULL),
      
      actionButton("run_gsea", "Run GSEA"),
      br(), br(),  # Adds space
      selectInput("logFC_col", "Select LogFC Column", choices = NULL),
      selectInput("pval_col", "Select P-Value Column", choices = NULL),   
      #selectInput("stat_col", "Select Stat Column", choices = NULL),
      actionButton("run_volcano", "Show Volcano Plot"),
      uiOutput("volcano_inputs"),
      hr(),  
      downloadButton("download_data", "Download dataset with inflammatome info"),
      downloadButton("download_plots", "Download plots")
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview",
                 DTOutput("table", width = "80%", height = "300px")
        ),
        tabPanel("GSEA Results", 
                 plotOutput("gsea_plot", height = "500px") %>% withSpinner(color="#0dc5c1"), 
                 DTOutput("gsea_table", width = "80%", height = "300px") 
        ),
        
        tabPanel("Volcano Plot", 
                 selectInput("volcano_set", "Choose inflammation gene set to highlight on plot:",
                             choices = c("Inflammatome (top 2000)" = "top2000", "Inflammation signature (top 100)" = "top100"),
                             selected = "top2000"),
                 plotOutput("volcano_plot", height = "500px"),
                 br(),
                 br(),
                 br(),
                 br(),
                 tags$div(style = "padding: 10px; border-top: 1px solid #ccc;",
                          tags$p(strong("Interpretation guide:"), style = "font-size: 16px;"),
                          tags$p("The vertical red line represents the median logFC of the genes from the selected inflammation gene set. When the line is shifted to the right compared to origin, it indicates some degree of inflammation is present in the dataset."),
                          tags$p("This plot can be interpreted or used in the following way: when many of the top differentially expressed genes are part of the inflammatome, one could consider filtering out those genes from further analysis or validation, in order to focus on more disease-specific genes and 
                                 not general inflammatory responses. For help in filtering inflammatome genes, you can download the dataset with the extra columns added for inflammatome info."),
                          style = "font-size: 14px; color: #333;"
                 )
        )
        #tabPanel("Volcano Plot", 
        #         plotOutput("volcano_plot", height = "500px"))
      )
    )),
    
    tags$head(
      tags$style(HTML("
    table.dataTable.compact tbody td {
      padding: 2px 6px !important;
      font-size: 12px !important
    }
  "))
 )
)


# Helper function to find the first matching column -------
select_best_match <- function(preferred, col_names, fallback = NULL) {
  match <- preferred[preferred %in% col_names]
  if (length(match) > 0) {
    return(match[1])  # Use the first match
  } else if (!is.null(fallback) && fallback %in% col_names) {
    return(fallback)  # Use fallback if available
  } else {
    return(col_names[1])  # Default to the first column
  }
}

# Server ---------
server <- function(input, output, session) {
  
  observeEvent(input$toggle_info, {
    shinyjs::toggle(id = "more_info_text", anim = TRUE, animType = "slide")
  })
  
  
  # Define example dataset choices
  example_datasets <- list(
    "Park et. al, 2022, IgA nephropathy, RNA-Seq" = "https://raw.githubusercontent.com/opalasca/inflammatome_package_sandbox/main/data/test_datasets/02_GSE175759_IgAN_ctl.tsv",
    "Bennike et. al, 2015, ulcerative colitis, proteomics" = "https://raw.githubusercontent.com/opalasca/inflammatome_package_sandbox/main/data/test_datasets/DE.res.UC.Andersen.raw.tsv"
  )
  
  # Reactive value to track whether example mode is active
  example_mode <- reactiveVal(FALSE)
  
  # Reactive value for storing dataset
  dataset <- reactiveVal(NULL)
  
  # Reactive value for storing selected sheet data from Excel
  selected_excel_data <- reactiveVal(NULL)
  
  # Reactive to handle the example mode
  observe({
    if (!is.null(input$example_file_choice) && input$example_file_choice != "") {
      example_mode(TRUE)  # Switch to example mode
      selected_url <- example_datasets[[input$example_file_choice]]
      showNotification(paste("Loading example dataset:", input$example_file_choice), type = "message", duration = 3)
      data <- read.delim(selected_url, header = TRUE, sep = "\t")
      dataset(data)
    } else {
      example_mode(FALSE)  # Switch to user-uploaded mode
    }
  })
  
  observe({
    req(input$file)
    
    # Reset selected dataset when switching file types
    dataset(NULL)
    selected_excel_data(NULL)
    
    # Check if file is an Excel file
    if (grepl("\\.xlsx$", input$file$name)) {
      example_mode(FALSE)  # Ensure we're not in example mode
      
      # Get the sheet names
      sheet_names <- readxl::excel_sheets(input$file$datapath)
      
      # Update sheet selection UI
      output$sheet_select_ui <- renderUI({
        selectInput("sheet_select", "Select Sheet", choices = sheet_names)
      })
      
    } else {
      ext <- tools::file_ext(input$file$datapath)
      
      if (ext == "csv") {
        data <- read.csv(input$file$datapath, header = TRUE)
      } else {
        data <- read.delim(input$file$datapath, header = TRUE, sep = "\t")
      }
      
      dataset(data)
      output$sheet_select_ui <- renderUI(NULL)
    }
  })
  
  # Reactive to load data from the selected Excel sheet
  observeEvent(input$sheet_select, {
    req(input$sheet_select, input$file)
    
    if (grepl("\\.xlsx$", input$file$name)) {
      # Read data from the selected sheet
      data <- readxl::read_excel(input$file$datapath, sheet = input$sheet_select)
      
      # Update both dataset() and selected_excel_data()
      selected_excel_data(data)
      dataset(data)  # Ensure dataset updates properly
    }
  })
  
  
  # Final dataset reactive that will combine TSV/TXT and Excel logic
  current_data <- reactive({
    if (example_mode()) {
      return(dataset())  # Return example dataset
    } else if (!is.null(dataset())) {
      return(dataset())  # Return TSV/TXT dataset
    } else if (!is.null(selected_excel_data())) {
      return(selected_excel_data())  # Return selected Excel sheet data
    }
    return(NULL)
  })
  
  # Observe dataset change to update UI components like column selections
  observeEvent(current_data(), {
    req(current_data())  # Ensure data exists
    data <- current_data()
    col_names <- colnames(data)
    
    # Update the select inputs for columns
    updateSelectInput(session, "id_col", choices = col_names, selected = col_names[1])
    updateSelectInput(session, "sort_col", choices = col_names, selected = select_best_match(c("stat", "t"), col_names, fallback = col_names[2]))
    updateSelectInput(session, "logFC_col", choices = col_names, selected = select_best_match(c("log2FoldChange", "logFC"), col_names))
    updateSelectInput(session, "pval_col", choices = col_names, selected = select_best_match(c("pvalue", "P.Value", "pval"), col_names))
    
    # Attempt to detect keytype
    selected_id_col <- col_names[1]
    id_vals <- as.character(data[[selected_id_col]])
    id_vals <- id_vals[!is.na(id_vals)]
    sample_vals <- head(id_vals, 100)
    
    detected_keytype <- "Symbol"  # Default fallback
    
    if (any(grepl("^ENSG[0-9]+", sample_vals))) {
      detected_keytype <- "Ensembl"
    } else if (any(grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", sample_vals))) {
      detected_keytype <- "Uniprot"
    } else if (all(grepl("^[0-9]+$", sample_vals))) {
      detected_keytype <- "Entrez"
    }
    
    updateSelectInput(session, "keytype", selected = detected_keytype)
  })
  
  # Processed Data
  processed_data <- reactive({
    req(current_data(), input$id_col, input$keytype)  # Ensure required inputs are available
    
    data <- current_data()
    id_col <- input$id_col
    keytype <- tolower(input$keytype)
    
    # Validate ID column
    if (!(id_col %in% colnames(data))) {
      showNotification("⚠️ Selected ID column does not exist!", type = "error", duration = 5)
      return(NULL)
    }
    
    # Try processing the data, handle keytype errors 
    result <- tryCatch({
      processed <- process_input_data(data, id = id_col, keytype = keytype)
      
      if (is.null(processed)) {
        showNotification("⚠️ Identifiers could not be mapped. Please check your keytype selection.", 
                         type = "warning", duration = 5)
        return(NULL)  # Allow re-selection instead of stopping
      }
      
      return(processed)
    }, error = function(e) {
      showNotification("⚠️ Invalid keytype selected. Please try again.", 
                       type = "error", duration = 5)
      return(NULL)  # Keep app running instead of stopping
    })
    
    return(result)
  })
 
  gene_sets <- reactive({
    file_path <- if (input$keytype == "Ensembl") {
      "data/gene_sets_Ensembl.tsv"
    } else {
      "data/gene_sets_Entrez.tsv"
    }
    
    # Read the gene sets file
    read_tsv(file_path, show_col_types = FALSE)  # show_col_types = FALSE prevents extra output
  })
  
  # GSEA Analysis
  gsea_results <- eventReactive(input$run_gsea, {
    req(processed_data(), input$sort_col)
    
    showNotification("Running GSEA analysis...", type = "message", duration = 5)
    
    data <- processed_data()
    id_col <- isolate(input$id_col)  # Ensure the input is only read when the button is clicked
    sort_col <- isolate(input$sort_col)
    keytype <- isolate(input$keytype)
    
    if (is.null(data)) {
      showNotification("⚠️ Processed data is NULL. Check keytype selection.", type = "error", duration = 5)
      return(NULL)
    }

    #gene_sets <- get_gene_sets(keytype)
    
    #result <- gsea_analysis(data, gene_sets(), sorting_value_col_name = sort_col, name = "data")
    result <- tryCatch({
      gsea_analysis(data, gene_sets(), sorting_value_col_name = sort_col, name = "data")
      }, error = function(e) {
      showNotification("⚠️ GSEA failed: Check keytype or sorting column.", type = "error", duration = 5)
      return(NULL)
    })
    
    showNotification("GSEA analysis complete!", type = "message", duration = 3)
    
    return(result)
  })
  
  # Show a notification and disable the button when running GSEA
  observeEvent(input$run_gsea, {
    shinyjs::disable("run_gsea")  # Disable button
    showNotification("Running GSEA analysis...", type = "message", duration = 5)
    
    # Trigger the GSEA analysis by calling gsea_results
    gsea_results()  
    
    # Show completion message and re-enable button after GSEA finishes
    showNotification("GSEA analysis complete!", type = "message", duration = 3)
    shinyjs::enable("run_gsea")  # Re-enable button after completion
  })
  
  # Enable the runGSEA button when keytype or other relevant inputs change
  observeEvent(input$keytype, {
    shinyjs::enable("run_gsea")
    data <- processed_data()  # Re-trigger data processing
  })
  
  observeEvent(input$sort_col, {
    shinyjs::enable("run_gsea")
  })
  
  observeEvent(input$id_col, {
    shinyjs::enable("run_gsea")
  })
  
  volcano_results <- eventReactive({
    input$run_volcano
    input$volcano_set  # make it reactive to dropdown
  }, {
    req(processed_data(), input$logFC_col, input$pval_col, input$volcano_set)
    
    data <- processed_data()
    
    if (input$volcano_set == "top2000") {
      plot_volcano_2000(
        data,
        keytype = input$keytype,
        logFC_col_name = input$logFC_col,
        pval_col_name = input$pval_col
      )
    } else {
      plot_volcano_100(
        data,
        keytype = input$keytype,
        logFC_col_name = input$logFC_col,
        pval_col_name = input$pval_col
      )
    }
  })
  
  
  # Volcano Plot
  #volcano_results <- eventReactive(input$run_volcano, {
  #  req(processed_data(), input$logFC_col, input$pval_col)
  #  data <- processed_data()
  #  plot_volcano(data, keytype=input$keytype, logFC_col_name = input$logFC_col, pval_col_name = input$pval_col)
  #})
  
  # Show a notification and disable the button when running volcano
  observeEvent(input$run_volcano, {
    shinyjs::disable("run_volcano")  # Disable button
    showNotification("Generating volcano plots...", type = "message", duration = 5)
    
    # Trigger the Volcano plot generation by calling gsea_results
    #volcano_results()  
    
    # Show completion message and re-enable button after GSEA finishes
    showNotification("Volcano plots complete!", type = "message", duration = 3)
    shinyjs::enable("run_volcano")  # Re-enable button after completion
  })
  
  
  # Render GSEA Plot
  output$gsea_plot <- renderPlot({
    req(gsea_results())
    plot_gsea(gsea_results(), "data")
  }, height = 500, width = 700)
  
  # Render GSEA Table
  output$gsea_table <- renderDT({
    req(gsea_results())
    datatable(gsea_results(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  
  # Render the dataset table
  output$table <- renderDT({
    req(current_data())  # Ensure data exists
    datatable(current_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Render processed data table
  output$processed_table <- renderDT({
    req(processed_data())  # Ensure processed data exists
    datatable(processed_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Render Volcano Plot
  output$volcano_plot <- renderPlot({
    req(volcano_results())
    volcano_results()
  }, height = 600, width = 700)
  
  # Download processed data  
  output$download_data <- downloadHandler(
    filename = function() { paste("dataset", Sys.Date(), ".csv", sep="") },
    content = function(file) {
      write.csv(processed_data(), file, row.names = FALSE)
    }
  )
  
  output$download_plots <- downloadHandler(
    filename = function() {
      paste("plots_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      # Create a subdirectory inside tempdir
      plot_dir <- file.path(tempdir(), "plots")
      dir.create(plot_dir, showWarnings = FALSE)
      
      # Define file paths
      gsea_plot_path <- file.path(plot_dir, "GSEA_plot.png")
      volcano_plot_path <- file.path(plot_dir, "Volcano_plot.png")
      
      # Save GSEA plot (if exists)
      if (!is.null(gsea_results())) {
        png(gsea_plot_path, width = 800, height = 600)
        plot_gsea(gsea_results(), "data")  
        dev.off()
      }
      
      # Save Volcano plot (if exists)
      if (!is.null(volcano_results())) {
        png(volcano_plot_path, width = 800, height = 450)
        print(volcano_results())  # Ensure ggplot gets printed
        dev.off()
      }
      
      # Get the list of saved plots
      plot_files <- list.files(plot_dir, full.names = TRUE)
      
      # Create a ZIP file only if plots exist
      if (length(plot_files) > 0) {
        zip(zipfile = file, files = normalizePath(plot_files))
      } else {
        showNotification("No plots available to download.", type = "error")
      }
    },
    contentType = "application/zip"
  )
  
  
}

# Run App
shinyApp(ui = ui, server = server)






























