library(shiny)
library(plotly)
library(Seurat)
library(ggplot2)
library(Matrix)

bin_mz <- function(data, mz_list, assay = "Spatial", slot = "counts", stored.in.metadata = FALSE){
  data_copy <- data

  if (stored.in.metadata){
    metadata_counts <- data_copy@meta.data[mz_list]
    if (length(colnames(metadata_counts)) < 2) {
      stop("One or more genes not found in the assay meta.data.")
    }
    binned.data <- Matrix::rowSums(metadata_counts)

  } else{
    assay_counts <- data_copy[[assay]][slot]
    selected_genes <- assay_counts[mz_list, , drop = FALSE]
    if (is.null(selected_genes)) {
      stop("One or more genes not found in the assay counts.")
    }
    binned.data <- Matrix::colSums(selected_genes)
  }

  return(binned.data)
}

binmetabolites <- function(data, mzs, assay = "Spatial", slot = "data", bin_name = "Binned_Metabolites"){

  binned_counts <- bin_mz(data, mz_list = mzs, assay = assay, slot = slot) # bins m/z masses

  data[[bin_name]]<- binned_counts #adds to metadata

  return(data)
}




bladder_HE <- readRDS("bladder_HE.RDS")
bladder_HE <- subset(bladder_HE, subset = ssc %in% c("2", "5", "6"))

obj <- bladder_HE
assay = "Spatial"
slot = "counts"
image = "slice1"

means <- Matrix::as.matrix(Matrix::rowMeans(obj[[assay]][slot]))
df <- data.frame(means)
df$x <- as.numeric(gsub("mz-", "", rownames(df)))
df$y <- df$means
rownames(df) <- NULL

# Shiny app UI
ui <- fluidPage(
  titlePanel("Interactive Spatial Intensity Plot"),
  sidebarLayout(
    sidebarPanel(
      numericInput("curve_center", "m/z Value:", value = 450, min = min(df$means), max = max(df$means), step = 1),
      radioButtons("bin_mode", "Binning Mode:",
                   choices = c("Absolute (m/z)" = "mz", "Relative (ppm)" = "ppm"),
                   selected = "mz"),
      fluidRow(
        column(6,
               sliderInput("curve_width_slider", "Bin Width:", min = 0, max = 100, value = 10, step = 0.1)
        ),
        column(6,
               numericInput("curve_width","", value = 10, min = 0, max = 1000)
        )
      ),
      numericInput("spot_size", "Plot Spot Size:", value = 10, min = 0, max = 100),  # Added spot size input

      fluidRow(
        column(6, numericInput("x_min", "X-Axis Minimum:", value = min(df$means), min = min(df$means), max = max(df$means), step = 1)),
        column(6, numericInput("x_max", "X-Axis Maximum:", value = max(df$means), min = min(df$means), max = max(df$means), step = 1))
      )
    ),
    mainPanel(
      plotlyOutput("plot"),
      verbatimTextOutput("selected_indices"),
      plotOutput("spatial_plot")
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$curve_width_slider, {
    updateNumericInput(session, "curve_width", value = input$curve_width_slider)
  })

  # Keep slider in sync with numeric input
  observeEvent(input$curve_width, {
    updateSliderInput(session, "curve_width_slider", value = input$curve_width)
  })


  output$plot <- renderPlotly({
    curve_center <- input$curve_center
    bin_mode <- input$bin_mode

    if (bin_mode == "ppm") {
      tolerance <- input$curve_width * 1e-6
      lower_bound <- curve_center - (curve_center * tolerance)
      upper_bound <- curve_center + (curve_center * tolerance)
      curve_width <- (upper_bound - lower_bound) / 2
    } else {
      curve_width <- input$curve_width
    }


    x_min <- input$x_min
    x_max <- input$x_max

    # Calculate normal distribution values based on 'means'
    curve_vals <- dnorm(df$means, mean = curve_center, sd = curve_width)

    # Scale curve to match max y value at curve_center
    max_y_at_center <- df$y[which.min(abs(df$means - curve_center))]

    # Generate new x-values for the curve
    new_x_vals <- seq(from = curve_center - curve_width, to = curve_center + curve_width, length.out = 100)
    new_curve_vals <- dnorm(new_x_vals, mean = curve_center, sd = curve_width)
    new_curve_vals <- (new_curve_vals - min(new_curve_vals)) / (max(new_curve_vals) - min(new_curve_vals))
    new_curve_vals <- new_curve_vals * max_y_at_center

    # Create histogram with curve overlay
    plot_ly() %>%
      add_bars(x = df$x, y = df$y, name = "Histogram") %>%
      add_lines(x = new_x_vals, y = new_curve_vals, name = "Curve", line = list(color = 'blue')) %>%
      layout(
        xaxis = list(title = "Means (m/z values)", range = c(x_min, x_max)),
        yaxis = list(title = "Frequency"),
        showlegend = TRUE
      )
  })

  selected_points <- reactive({
    curve_center <- input$curve_center
    curve_width <- input$curve_width
    bin_mode <- input$bin_mode

    if (bin_mode == "ppm") {
      tolerance <- curve_width * 1e-6
      lower_bound <- curve_center - (curve_center * tolerance)
      upper_bound <- curve_center + (curve_center * tolerance)
      curve_width <- (upper_bound - lower_bound) / 2
    } else {
      lower_bound <- curve_center - curve_width
      upper_bound <- curve_center + curve_width
    }

    df$x[df$x >= lower_bound & df$x <= upper_bound]
  })

  output$selected_indices <- renderText({
    indices <- selected_points()

    if (length(indices) > 0) {
      paste("Selected feature:", paste(paste0("mz-", indices), collapse = ", "))
    } else {
      "No rows under the curve."
    }
  })

  output$spatial_plot <- renderPlot({
    indices <- selected_points()
    spot_size <- input$spot_size

    if (length(indices) > 0) {

      curve_center <- input$curve_center
      curve_width <- input$curve_width
      bin_mode <- input$bin_mode

      if (bin_mode == "ppm") {
        tolerance <- curve_width * 1e-6
        lower_bound <- curve_center - (curve_center * tolerance)
        upper_bound <- curve_center + (curve_center * tolerance)
        curve_width <- (upper_bound - lower_bound) / 2
      }

      features_to_plot <- indices  # Example gene names based on index
      features_to_plot <- paste0("mz-",features_to_plot)
      feature_name <- sprintf("mz-%d \u00B1 %.1f", curve_center, curve_width)
      binned_obj <- binmetabolites(obj, mzs = features_to_plot, assay = "Spatial", slot = "counts", bin_name = feature_name)
      Seurat::SpatialFeaturePlot(binned_obj, images = image, features = feature_name,  pt.size.factor = spot_size) &theme_void()
    } else {
      plot(1, type = "n", xlab = "", ylab = "", main = "No selected features")
      text(1, 1, "No data available")
    }
  })



}

shinyApp(ui, server)
