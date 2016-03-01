library(shiny)

shinyUI(fluidPage(
    titlePanel("Interactive Visualization of Cytofkit Results"),
    
    fluidRow(
        column(3,
               
               fileInput('cytofkitObj', strong('Cytofkit RData:'), multiple = FALSE, 
                         accept = c('text/RData', '.RData')),
               actionButton("goButton", "Submit!"),
               
               hr(),
               uiOutput("sampleSelect"),
               
               hr(),
               h5("If Add Cluster Labels"),
               checkboxInput("addLabel", label = "Add Cluster Label", value = TRUE),
               
               hr(),
               h4("Summary:"),
               h5("Expression Data:"),
               textOutput("summaryText1"),
               h5("Cluster Method(s):"),
               textOutput("summaryText2"),
               h5("Visualization Method(s):"),
               textOutput("summaryText3"),
               h5("Progression Method(s):"),
               textOutput("summaryText4"),
               
               hr(),
               div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
               div(style = "margin-top: 10px; ", 
                   HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
        ),
        column(9,
               tabsetPanel(type = "pills",
                           
                           tabPanel("Scatter Plot", fluidPage(
                               hr(),
                               fluidRow(
                                   column(2,
                                          uiOutput("x_choose")
                                   ),
                                   column(2, 
                                          uiOutput("y_choose")
                                   ),
                                   column(4, 
                                          uiOutput("z_choose")
                                   ),
                                   column(2,
                                          numericInput("pointSize", "Point Size:", value = 3)
                                   ),
                                   column(2, 
                                          numericInput("labelSize", "Label Size:", value = 12)
                                   )
                               ),
                               
                               hr(),
                               plotOutput("xyzPlot", width = "80%")
                           )),
                           
                           tabPanel("Heatmap", fluidPage(
                               hr(),
                               fluidRow(
                                   column(4, 
                                          uiOutput("cMethod_choose")
                                   ),
                                   column(4,
                                          selectInput('heatmapMethod', strong('Heatmap Type:'), 
                                                      choices = c("mean", "median", "percentage"), 
                                                      selected = "mean", width = "100%")
                                   ),
                                   column(2,
                                          numericInput("rLabelSize", "Row Label Size:", value = 1, step = 0.5)
                                   ),
                                   column(2, 
                                          numericInput("cLabelSize", "Col Label Size:", value = 1, step = 0.5)
                                   )
                               ),
                               
                               hr(),
                               plotOutput("heatmapPlot", width = "100%")
                               
                           )),
                           
                           tabPanel("Progression", fluidPage(
                               hr(),
                               uiOutput("pMethod_choose"),
                               hr(),
                               plotOutput("progressionPlot", width = "100%")
                           )) 
            )
        )
    )
))