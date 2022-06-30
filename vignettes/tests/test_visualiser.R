library(shiny)

ui <- basicPage(sliderInput("range", "Age:",min = 0, max = 100, value = c(0,100)),textOutput("SliderText"))
server <- shinyServer(function(input, output, session){
  my_range <- reactive({
    cbind(input$range[1],input$range[2])
  })
  output$SliderText <- renderText({my_range()})
})
shinyApp(ui = ui, server = server)