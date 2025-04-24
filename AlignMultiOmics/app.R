library(shiny)
library(shinyjs)
library(plotly)
library(Seurat)
library(ggplot2)
library(Matrix)
library(zeallot)


################################################################################
sm.data <- readRDS("sm.RDS")
st.data <- readRDS("st.RDS")

cat("loaded data")
rotate <- function (
    angle,
    center.cur
) {
  alpha <- 2*pi*(angle/360)
  #center.cur <- c(center.cur, 0)
  #points(center.cur[1], center.cur[2], col = "red")
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.transf(center.cur[1], center.cur[2], alpha)%*%tr
  return(tr)
}


#' Creates a transformation matrix that translates an object
#' in 2D
#'
#' @param translate.x,translate.y translation of x, y coordinates

translate <- function (
    translate.x,
    translate.y
) {
  tr <- rigid.transl(translate.x, translate.y)
  return(tr)
}


#' Creates a transformation matrix that mirrors an object
#' in 2D along either the x axis or y axis around its
#' center of mass
#'
#' @param mirror.x,mirror.y Logical specifying whether or not an
#' object should be reflected
#' @param center.cur Coordinates of the current center of mass
#'

mirror <- function (
    mirror.x = FALSE,
    mirror.y = FALSE,
    center.cur
) {
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.refl(mirror.x, mirror.y)%*%tr
  tr <- rigid.transl(center.cur[1], center.cur[2])%*%tr
  return(tr)
}


#' Stretch along angle
#'
#' Creates a transformation matrix that stretches an object
#' along a specific axis
#'
#' @param r stretching factor
#' @param alpha angle
#' @param center.cur Coordinates of the current center of mass
#'

stretch <- function(r, alpha, center.cur) {
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.rot(alpha, forward = TRUE)%*%tr
  tr <- rigid.stretch(r)%*%tr
  tr <- rigid.rot(alpha, forward = FALSE)%*%tr
  tr <- rigid.transl(center.cur[1], center.cur[2])%*%tr
  return(tr)
}


#' Creates a transformation matrix for rotation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#'
#' @param alpha rotation angle
#' @param forward should the rotation be done in forward direction?
#'

rigid.rot <- function (
    alpha = 0,
    forward = TRUE
) {
  alpha <- 2*pi*(alpha/360)
  tr <- matrix(c(cos(alpha), ifelse(forward, -sin(alpha), sin(alpha)), 0, ifelse(forward, sin(alpha), -sin(alpha)), cos(alpha), 0, 0, 0, 1), nrow = 3)
  return(tr)
}


#' Creates a transformation matrix for rotation and translation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#' followed by a translation with an offset of (h, k). Points are assumed to be
#' centered at (0, 0).
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#' @param alpha rotation angle
#'

rigid.transf <- function (
    h = 0,
    k = 0,
    alpha = 0
) {
  tr <- matrix(c(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, h, k, 1), nrow = 3)
  return(tr)
}

#' Creates a transformation matrix for translation with an offset of (h, k)
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#'

rigid.transl <- function (
    h = 0,
    k = 0
) {
  tr <-  matrix(c(1, 0, 0, 0, 1, 0, h, k, 1), nrow = 3)
  return(tr)
}

#' Creates a transformation matrix for reflection
#'
#' Creates a transformation matrix for reflection where mirror.x will reflect the
#' points along the x axis and mirror.y will reflect thepoints along the y axis.
#' Points are assumed to be centered at (0, 0)
#'
#' @param mirror.x,mirror.y Logical: mirrors x or y axis if set to TRUE

rigid.refl <- function (
    mirror.x,
    mirror.y
) {
  tr <- diag(c(1, 1, 1))
  if (mirror.x) {
    tr[1, 1] <- - tr[1, 1]
  }
  if (mirror.y) {
    tr[2, 2] <- - tr[2, 2]
  }
  return(tr)
}

#' Creates a transformation matrix for stretching
#'
#' Creates a transformation matrix for stretching by a factor of r
#' along the x axis.
#'
#' @param r stretching factor

rigid.stretch <- function (
    r
) {
  tr <- matrix(c(r, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
}


#' Combines rigid tranformation matrices
#'
#' Combines rigid tranformation matrices in the following order:
#' translation of points to origin (0, 0) -> reflection of points
#' -> rotation by alpha degrees and translation of points to new center
#'
#' @param center.cur (x, y) image pixel coordinates specifying the current center of the tissue (stored in slot "tools" as "centers")
#' @param center.new (x, y) image pixel coordinates specifying the new center (image center)
#' @param alpha Rotation angle
#'
#' @inheritParams rigid.transf
#' @inheritParams rigid.transl
#' @inheritParams rigid.refl
#'
#' @examples
#' \dontrun{
#' library(imager)
#' library(tidyverse)
#' im <- load.image("https://upload.wikimedia.org/wikipedia/commons/thumb/f/fd/Aster_Tataricus.JPG/1024px-Aster_Tataricus.JPG")
#' d <- sRGBtoLab(im) %>% as.data.frame(wide="c")%>%
#'   dplyr::select(-x,-y)
#'
#' km <- kmeans(d, 2)
#'
#' # Run a segmentation to extract flower
#' seg <- as.cimg(abs(km$cluster - 2), dim = c(dim(im)[1:2], 1, 1))
#' plot(seg); highlight(seg == 1)
#'
#' # Detect edges
#' dx <- imgradient(seg, "x")
#' dy <- imgradient(seg, "y")
#' grad.mag <- sqrt(dx^2 + dy^2)
#' plot(grad.mag)
#'
#' # Extract points at edges
#' edges.px <- which(grad.mag > max(grad.mag[, , 1, 1])/2, arr.ind = TRUE)
#' points(edges.px, col = "green", cex = 0.1)
#'
#' # Apply transformations to point set
#' tr1 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(1200, 1200), alpha = 90)
#' tr2 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(500, 1200), mirror.x = T, alpha = 30)
#' tr3 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(1200, 500), mirror.y = T, alpha = 270)
#' plot(edges.px, xlim = c(0, 1700), ylim = c(0, 1700), cex = 0.1)
#' points(t(tr1%*%t(edges.px[, 1:3])), cex = 0.1, col = "red")
#' points(t(tr2%*%t(edges.px[, 1:3])), cex = 0.1, col = "yellow")
#' points(t(tr3%*%t(edges.px[, 1:3])), cex = 0.1, col = "blue")
#' }
#'
#' @export

combine.tr <- function (
    center.cur,
    center.new,
    alpha,
    mirror.x = FALSE,
    mirror.y = FALSE
) {

  alpha <- 2*pi*(alpha/360)
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])

  # reflect
  tr <- rigid.refl(mirror.x, mirror.y)%*%tr

  # rotate and translate to new center
  tr <- rigid.transf(center.new[1], center.new[2], alpha)%*%tr
}



generate.map.affine <- function (
    tr,
    forward = FALSE
) {

  if (forward) {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      xy <- t(solve(tr)%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  } else {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      xy <- t(tr%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  }
  return(map.affine)
}


## ST cooridnates
df <- st.data[c("x", "y")]

## SM Coordinates
df2 <- sm.data[c("x", "y")]
df2$x <- df2$x * 20 / (10/3)
df2$y <- df2$y * 20 / (10/3)
df2 <- df2[c("x", "y")]

# Calculate scatter for plotting
df$pixel_x <- df$x
df$pixel_y <- df$y

sc1 <- df[c("x", "y")]
rownames(sc1) <- NULL
coords1 <- df[c("pixel_x", "pixel_y")]

df2$pixel_x <- df2$x
df2$pixel_y <- df2$y

sc2 <- df2[c("x", "y")]
rownames(sc2) <- NULL
coords2<- df2[c("pixel_x", "pixel_y")]


counter <- NULL
coords.ls <- NULL
scatters <- list("1" = list("scatter" = sc1, "coords" = coords1),
                 "2" = list("scatter" = sc2, "coords" = coords2))
rm(coords2)
rm(sc2)
rm(sc1)
rm(coords1)
rm(df2)
rm(df)


tr.matrices <- lapply(list(diag(c(1, 1, 1)), diag(c(1, 1, 1))), function(x) diag(c(1, 1, 1)))

image.dims <- list("1" = c(594,600), "2" = c(594,600))

sm.data$x <- NULL
sm.data$y <- NULL

st.data$x <- NULL
st.data$y <- NULL
################################################################################


cat("running")


ui <- fluidPage(
  useShinyjs(),
  fluidRow(
    column(5,
           shiny::hr(),
           actionButton(inputId = "info", label = "Instructions"),
           shiny::hr(),
           # Combine sliders for shifting and stretching in fewer rows
           fluidRow(
             column(6, sliderInput("angle", "Rotation angle", value = 0, min = -120, max = 120, step = 0.1)),
             column(6, sliderInput("shift_x", "Move along x axis", value = 0, min = -450, max = 450, step = 1)),
             column(6, sliderInput("shift_y", "Move along y axis", value = 0, min = -450, max = 450, step = 1)),
             column(6, sliderInput("stretch_angle1", "Stretch angle 1", value = 0, min = -180, max = 180, step = 0.1)),
             column(6, sliderInput("stretch_factor1", "Stretch factor 1", value = 1, min = 0.1, max = 2, step = 0.01)),
             column(6, sliderInput("stretch_angle2", "Stretch angle 2", value = 0, min = -180, max = 180, step = 0.1)),
             column(6, sliderInput("stretch_factor2", "Stretch factor 2", value = 1, min = 0.1, max = 2, step = 0.01))
           ),
           # Group numeric inputs
           fluidRow(
             column(4, numericInput("size_spot", "SM spot size", value = 0.5, min = 0, max = 5, step = 0.1)),
             column(4, numericInput("size_target", "ST point size", value = 0.3, min = 0, max = 5, step = 0.05)),
             column(4, selectInput("spot_shape", "Spot shape", choices = c("circle", "square")))
           ),
           # Group feature selectors for better layout
           fluidRow(
             column(4, selectInput("sm_plot", "SM plot feature", choices = colnames(sm.data))),
             column(4, selectInput("st_plot", "ST plot feature", choices = colnames(st.data)))
           ),
           # Group checkboxes together
           fluidRow(
             column(4, checkboxInput("show_ST_img", "Show image", value = TRUE)),
             column(4, checkboxInput("show_ST_spots", "Show ST spots", value = FALSE)),
             column(4, checkboxInput("show_SM_spots", "Show SM spots", value = TRUE)),
             column(4, checkboxInput("flip_x", "Mirror along x axis", value = FALSE)),
             column(4, checkboxInput("flip_y", "Mirror along y axis", value = FALSE))
           ),
           actionButton("myBtn", "Return aligned data")
    ),
    column(7, plotOutput("scatter", width = "100%", height = "100%"))
  )
)

server <- function(input, output) {

  rotation_angle <- reactive({
    rot_angle <- input$angle
    return(rot_angle)
  })

  translation_xy <- reactive({
    trxy <- c(input$shift_x, input$shift_y)
    return(trxy)
  })

  mirror_xy <- reactive({
    mirrxy <- c(input$flip_x, input$flip_y)
    return(mirrxy)
  })

  stretch_angle1 <- reactive({
    str_angle1 <- input$stretch_angle1
    return(str_angle1)
  })

  stretch_factor1 <- reactive({
    str_factor1 <- input$stretch_factor1
    return(str_factor1)
  })

  stretch_angle2 <- reactive({
    str_angle2 <- input$stretch_angle2
    return(str_angle2)
  })

  stretch_factor2 <- reactive({
    str_factor2 <- input$stretch_factor2
    return(str_factor2)
  })


  pt_size <- reactive({
    input$size_spot
  })

  st_feature_plot <- reactive({
    input$st_plot
  })

  sm_feature_plot <- reactive({
    input$sm_plot
  })

  pt_shape <- reactive({
    if (input$spot_shape == "square"){
      return(15)
    } else{
      return(19)
    }
  })

  pt_size_target <- reactive({
    input$size_target
  })


  pt_st_points <- reactive({
    input$show_ST_spots
  })

  pt_sm_points <- reactive({
    input$show_SM_spots
  })

  pt_st_img <- reactive({
    input$show_ST_img
  })



  coords_list <- reactive({

    # Obtain point set and spot pixel coordinates
    ls <- scatter.coords()
    scatter.t <- ls[[1]]; coords.t <- ls[[2]]

    # Set transformation parameters
    xt.yt <- translation_xy()
    xy.alpha <- rotation_angle()
    mirrxy <-  mirror_xy()
    str.alpha1 <- stretch_angle1()
    str.factor1 <- stretch_factor1()
    str.alpha2 <- stretch_angle2()
    str.factor2 <- stretch_factor2()

    # Apply reflections
    center <- apply(scatter.t, 2, mean)
    tr.mirror <- mirror(mirror.x = mirrxy[1], mirror.y = mirrxy[2], center.cur = center)

    # Apply rotation
    tr.rotate <- rotate(angle = -xy.alpha, center.cur = center)

    # Apply translation
    tr.translate <- translate(translate.x = xt.yt[1], translate.y = -xt.yt[2])

    # Apply stretch
    tr.stretch1 <- stretch(r = str.factor1, alpha = -str.alpha1, center.cur = center)
    tr.stretch2 <- stretch(r = str.factor2, alpha = -(str.alpha2 + 90), center.cur = center)

    # Combine transformations
    tr <- tr.stretch2%*%tr.stretch1%*%tr.translate%*%tr.rotate%*%tr.mirror


    # Apply transformations
    scatter.t <- t(tr%*%rbind(t(scatter.t), 1))[, 1:2]
    coords.t <- t(tr%*%rbind(t(coords.t), 1))[, 1:2]

    return(list(scatter = scatter.t, coords = coords.t, tr = tr, xylimits = image.dims[[2]]))
  })

  output$scatter <- renderPlot({

    coords.ls <<- coords_list()
    c(scatter.t, coords.t, tr, xylimit) %<-% coords.ls

    d <- round((sqrt(xylimit[1]^2 + xylimit[2]^2) - xylimit[2])/2)

    center <- apply(coords.t[, 1:2], 2, mean)

    arrows.1 <- function(x0, y0, length.ar, angle.ar, ...){

      angle.ar <- 2*pi*(-angle.ar/360)
      ab <- cos(angle.ar) * length.ar
      bc <- sign(sin(angle.ar)) * sqrt(length.ar^2 - ab^2)

      x1 <- x0 + ab
      y1 <- y0 + bc

      arrows(x0, y0, x1, y1, ...)
    }



    cont_pal  <- RColorBrewer::brewer.pal("Reds", n = 9)
    cat_pal <- c("black", RColorBrewer::brewer.pal("Paired", n = 10))



    if (is.numeric(sm.data[[sm_feature_plot()]])){
      sm_cols <- cont_pal[as.numeric(cut(sm.data[[sm_feature_plot()]],breaks = 9))]
    } else {
      sm_cols <- cat_pal[as.factor(sm.data[[sm_feature_plot()]])]
    }

    if (is.numeric(st.data[[st_feature_plot()]])){
      st_cols <- cont_pal[as.numeric(cut(st.data[[st_feature_plot()]],breaks = 9))]
    } else {
      st_cols <- cat_pal[as.factor(st.data[[st_feature_plot()]])]
    }


    if (pt_st_img()){
      plot(readRDS("img.RDS"))

    } else{

      plot(NULL, NULL, col = "white", xlim = c(0, 594), ylim = c(0, 600), xaxt = 'n', yaxt = 'n', ann = FALSE, bty = "n")
    }

    if (pt_st_points()){
      points(scatters[[1]]$scatter[, 1], scatters[[1]]$scatter[, 2], col = st_cols, pch = as.numeric(pt_shape()), cex = pt_size_target())
    }

    if (pt_sm_points()){
      points(coords.t[, 1], coords.t[, 2], col = sm_cols, pch = as.numeric(pt_shape()), cex = pt_size())
      arrows.1(x0 = center[1], y0 = center[2], angle.ar = stretch_angle1(), length.ar = 100*stretch_factor1(), lwd = 4, col = "blue")
      arrows.1(x0 = center[1], y0 = center[2], angle.ar = 90 + stretch_angle2(), length.ar = 100*stretch_factor2(), lwd = 4, col = "red")
    }

  }, height = 500, width = 500)

  scatter.coords <- eventReactive(2, {
    reset("angle"); reset("shift_x"); reset("shift_y"); reset("flip_x"); reset("flip_y"); reset("stretch_factor1"); reset("stretch_factor2"); reset("stretch_angle1"); reset("stretch_angle2")
    if (!is.null(counter)) {
      scatters[[counter]] <<- coords.ls[c(1, 2)]
      if (!is.null(tr.matrices[[counter]])) {
        tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
      } else {
        tr.matrices[[counter]] <<- coords.ls[[3]]
      }
    }
    scatter <- scatters[[as.numeric(2)]]$scatter
    coords <- scatters[[as.numeric(2)]]$coords
    counter <<- as.numeric(2)
    return(list(scatter, coords))
  })

  observe({
    if(input$myBtn > 0){
      if (!is.null(counter)) {
        scatters[[counter]] <<- coords.ls[c(1, 2)]
        if (!is.null(tr.matrices[[counter]])) {
          tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
          cat("Sample:", counter, "\n",  tr.matrices[[counter]][1, ], "\n", tr.matrices[[counter]][2, ], "\n", tr.matrices[[counter]][3, ], "\n\n")
        } else {
          tr.matrices[[counter]] <<- coords.ls[[3]]
        }
      }
      stopApp(tr.matrices)
    }
  })

  observeEvent(input$info, {
    showModal(modalDialog(
      title = "Instructions",
      HTML(
        "The alignment interface is modified from [STUtility](https://github.com/jbergenstrahle/STUtility/tree/master) to provide an interface for manually aligning SpaMTP (Seurat) Objects <br>",
        "This interface can be used for: <br>",
        "1. Aligning SM data to ST data coordinates.<br><br>",
        "2. Aligning SM data to a H&E image <br>",
        "<br>",
        "How to use: <br>",
        "- Adjust the coordinates of the SM data by changing the rotation, x-axis or y-axis position to match the provided reference dataset. <br>",
        "- Stretch the SM coordinates in either the direction matching the blue and/or red arrow. The angle of the arrows can also be change to match the required stretching direction. <br>",
        "- Once them SM data is aligned select the 'Return Aligned Data' button to generate a SpaMTP Seurat object with the adjusted coordinates. <br>",
        "<br>",
        "Note: If using in 'AlignSpatialOmics' mode then a SpaMTP object will be returned containing only the SM data with the ajusted coordinates. <br>",
        "For mapping SM data to corresponding ST spots, please run 'MapSpatialOmics()' with the original ST object and now updated SM object. <br>"
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
}

# Returned transformation matrices
shinyApp(ui, server)


