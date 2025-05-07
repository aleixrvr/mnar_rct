library(shiny)
library(DiagrammeR)

ui <- fluidPage(

  tags$style(HTML(
    "
    .label-left .form-group {
      display: flex;              /* Use flexbox for positioning children */
      flex-direction: row;        /* Place children on a row (default) */
      width: 100%;                /* Set width for container */
      max-width: 400px;
    }
    .label-left label {
      margin-right: 2rem;         /* Add spacing between label and slider */
      align-self: center;         /* Vertical align in center of row */
      text-align: right;
      flex-basis: 100px;          /* Target width for label */
    }
    .label-left .irs { flex-basis: 300px;}  /* Target width for slider */
     "
  )),
  
  
# Application title

titlePanel("Missings"),

    sidebarLayout(
    sidebarPanel(
    div(class = "label-left",
     sliderInput("t_o", "T -> O", min = -1, max = 1, value = 0, step=0.01),
     sliderInput("o_a", "O ->  A", min = -1, max = 1, value = 0, step=0.01),
     sliderInput("x_o", "X ->  O", min = -1, max = 1, value = 0, step=0.01),
     sliderInput("x_a", "X ->  A", min = -1, max = 1, value = 0, step=0.01)
      )) ,
    
    
    mainPanel(
      
      fluidRow(
        column(6, 
               h3("Ground Truth: (O)"),
               textOutput("ATE1"),
               textOutput("EF11"),
               textOutput("EF12"),
        ),
        column(6,       
               h3("Raw estimates: (O*)"),
               textOutput("ATE2"),
               textOutput("EF21"),
               textOutput("EF22"), )
      ),
           h3("DAG: "),
           grVizOutput('DAG', width = "70%", height = "250px") 
        )
    )
)

server <- function(input, output) {

T0A1 <- reactive({  ((0.5+input$o_a/2) + input$x_a/2)*((0.5-input$t_o/2) + input$x_o/2) +
    ((0.5-input$o_a/2) + input$x_a/2)*(1- ((0.5-input$t_o/2) + input$x_o/2)) + 
    (0.5+input$o_a/2 - input$x_a/2)*(0.5-input$t_o/2 -input$x_o/2) + 
    (0.5-input$o_a/2 - input$x_a/2)*(1-(0.5-input$t_o/2 -input$x_o/2));
})
  
T1A1 <- reactive({ ((0.5+input$o_a/2) + input$x_a/2)*((0.5+input$t_o/2) + input$x_o/2) +
    ((0.5-input$o_a/2) + input$x_a/2)*(1- ((0.5+input$t_o/2) + input$x_o/2)) + 
    (0.5+input$o_a/2 - input$x_a/2)*(0.5+input$t_o/2 -input$x_o/2) + 
    (0.5-input$o_a/2 - input$x_a/2)*(1-(0.5+input$t_o/2 -input$x_o/2));
})

efecte0 <- reactive({ (0.5+input$o_a/2+input$x_a/2)*(0.5-input$t_o/2+input$x_o/2) / T0A1() + (0.5+input$o_a/2-input$x_a/2)*(0.5-input$t_o/2-input$x_o/2) / T0A1() ; })
efecte1 <- reactive({(0.5+input$o_a/2+input$x_a/2)*(0.5+input$t_o/2+input$x_o/2) / T1A1() + (0.5+input$o_a/2-input$x_a/2)*(0.5+input$t_o/2-input$x_o/2) / T1A1() ; })
  
  output$ATE1 <- renderText({ paste("Ground Truth ATE:",round(input$t_o,2))})
  output$EF11 <- renderText({ paste("G1:",round(0.5+input$t_o/2,2))})
  output$EF12 <- renderText({ paste("G2:",round(0.5-input$t_o/2,2))})

  output$ATE2 <- renderText({ paste("Biased ATE :", round(efecte1() - efecte0(),3))})
  output$EF21 <- renderText({ paste("G1:",round(efecte1(),2))})
  output$EF22 <- renderText({ paste("G2:",round(efecte0(),2))})
  
    
  output$DAG <- renderGrViz({ 
    DAG0 <- "digraph {
  A;O;X;T;Y[color = blue,label='O*'];
  {rank=same; T; X}
  {rank=same; O; A}"
    a1 <- paste('X -> O[label="',input$x_o,'"]')
    a2 <- paste('A -> Y[label="',"*",'"]')
    a3 <- paste('O -> A[label="',input$o_a,'"]')
    a4 <- paste('O -> Y[label="',"1",'"]')
    a5 <- paste('T -> O[label="',input$t_o,'"]')
    a6 <- paste('X -> A[label="',input$x_a,'"]')
    DAG1 <- "}";
    grViz(paste(DAG0,a1,a2,a3,a4,a5,a6, DAG1, sep=";"))  })

  }

shinyApp(ui = ui, server = server)






