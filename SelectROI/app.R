library(shiny)
library(plotly)
library(Seurat)
library(ggplot2)

seurat_obj <- read.csv("./spotted_roi_df.csv", row.names = 1)


coords <- seurat_obj[c("x","y")]


meta_cols <- colnames(seurat_obj)


ui = fluidPage(
  titlePanel("Select Regions of Interest"),
  sidebarLayout(
    sidebarPanel(
      selectInput("meta_col", "Select Metadata Column to Plot", choices = meta_cols),
      radioButtons("plot_type", "Plot Type",
                   choices = c("Feature (Continuous)" = "feature", "Categorical" = "categorical")),
      sliderInput("pt_size", "Spot Size", min = 1, max = 20, value = 5),
      textInput("roi_name", "Name for ROI", value = "ROI_1"),
      actionButton("save_roi", "Save ROI"),
      actionButton("reset_roi", "Reset ROI Selection"),
      actionButton("done", "Finish & Return Object"),
      verbatimTextOutput("status")
    ),
    mainPanel(
      plotlyOutput("spatial_plot", height = "600px")
    )
  )
)

server = function(input, output, session) {
  rv <- reactiveVal(seurat_obj)
  roi_mask <- reactiveVal(rep(0, nrow(coords)))
  
  output$spatial_plot <- renderPlotly({
    meta_col <- input$meta_col
    plot_type <- input$plot_type
    meta_data <- rv()[[meta_col]]
    spot_size <- input$pt_size
    
    plot_ly() %>%
      add_trace(
        x = coords$x,
        y = coords$y,
        type = "scattergl",
        mode = "markers",
        marker = if (plot_type == "feature") {
          list(color = meta_data, colorscale = "Viridis", showscale = TRUE, size = spot_size)
        } else {
          list(color = as.factor(meta_data), showscale = FALSE, size = spot_size)
        },
        text = rownames(coords),
        hoverinfo = "text"
      ) %>%
      layout(
        title = paste("Spatial Plot:", meta_col),
        xaxis = list(title = "X Coordinate"),
        yaxis = list(title = "Y Coordinate"),
        dragmode = "lasso"
      )
  })
  
  observeEvent(event_data("plotly_selected"), {
    sel_data <- event_data("plotly_selected")
    if (!is.null(sel_data)) {
      sel_points <- sel_data$pointNumber + 1
      sel_coords <- coords[sel_points, ]
      if (nrow(sel_coords) >= 3) {
        sel_coords <- sel_coords[chull(sel_coords[, c("x", "y")]), ]
        poly_coords <- as.matrix(rbind(sel_coords[, c("x", "y")], sel_coords[1, c("x", "y")]))
        poly <- st_polygon(list(poly_coords))
        poly_sf <- st_sfc(poly)
        poly_sf <- st_make_valid(poly_sf)
        
        points_sf <- st_as_sf(coords, coords = c("x", "y"))
        inside <- st_within(points_sf, poly_sf, sparse = FALSE)[, 1]
        
        new_mask <- roi_mask()
        new_mask[inside] <- 1
        roi_mask(new_mask)
        output$status <- renderText("Points selected for ROI.")
      } else {
        showNotification("Select at least 3 points", type = "warning")
      }
    }
  })
  
  observeEvent(input$save_roi, {
    name <- input$roi_name
    if (name == "") {
      showNotification("Please enter a name for the ROI.", type = "error")
      return()
    }
    
    seurat <- rv()
    seurat[[name]] <- roi_mask()
    rv(seurat)
    roi_mask(rep(0, nrow(coords)))
    output$status <- renderText(paste0("Saved ROI to metadata column: ", name))
  })
  
  observeEvent(input$reset_roi, {
    roi_mask(rep(0, nrow(coords)))
    output$status <- renderText("ROI selection reset.")
  })
  
  observeEvent(input$done, {
    stopApp(rv())
  })
}

shinyApp(ui, server)

