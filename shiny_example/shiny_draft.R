library(shiny)
setwd(here::here())

formatLambda <- function(lambda) {
  if (lambda == 0) {
    return("0.00")
  } else if ( (lambda * 1000) %% 10 == 0) {
    return(sprintf("%.2f", lambda))
  } else {
    return(sprintf("%.3f", lambda))
  }
}

ui <- fluidPage(
  titlePanel("L1 logistic result"),
  sidebarLayout(
    sidebarPanel(
      selectInput("eta", "select Eta:",
                  choices = c("1.1" = "1.1", "1.2" = "1.2")),
      selectInput("category", "select category:",
                  choices = c("alpha", "acinar", "beta", "delta", "ductal", "gamma", "class_combined")),
      sliderInput("lambda", "select Lambda:",
                  min = 0.02, max = 0.04, value = 0.03, step = 0.01)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("category predict result ", imageOutput("categoryImg")),
        tabPanel("Non zero parameters", imageOutput("combinedPlotImg")),
        tabPanel("confusion matrix", tableOutput("confMat"))
      )
    )
  )
  #uiOutput("dynamicImg")
)

server <- function(input, output) {

  output$categoryImg <- renderImage({
    eta <- input$eta
    category <- input$category
    dir <- paste0(here::here(),"/shiny_example/eta_", eta,"//")
    filePath <- file.path(dir, paste0(category, ".png"))

    if(eta == "1.2")
    {
      stop(paste(" Sorry Eta equals 1.2 doesn`t have category predict result till now "))
    }

    if (!file.exists(filePath)) {
      stop(paste("The file", filePath, "does not exist."))
    }

    list(src = filePath, contentType = "image/png", alt = "Category image.",
         width = "100%", height = "auto")
  }, deleteFile = FALSE)

  output$combinedPlotImg <- renderImage({
    eta <- input$eta
    lambda <- formatLambda(input$lambda)
    dir <- paste0(here::here(),"/shiny_example/eta_", eta, "/result_", lambda)
    filePath <- file.path(dir, "combined_plot.png")

    if (!file.exists(filePath)) {
      stop(paste("The file", filePath, "does not exist."))
    }

    list(src = filePath, contentType = "image/png", alt = "Combined plot image.",
         width = "100%", height = "auto")
  }, deleteFile = FALSE)

  output$confMat <- renderTable({
    eta <- input$eta
    lambda <- formatLambda(input$lambda)
    dir <- paste0(here::here(),"/shiny_example/eta_", eta, "/result_", lambda)
    filePath <- file.path(dir, "confusion_matrix.csv")

    if (!file.exists(filePath)) {
      stop(paste("The file", filePath, "does not exist."))
    }

    read.csv(filePath)
  })
}

shinyApp(ui = ui, server = server)
