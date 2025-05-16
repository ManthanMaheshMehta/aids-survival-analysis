library(shiny)
library(survival)
library(ggplot2)


aids_data <- read.csv("AIDS_Clinical_Trials_Study_175.csv")


cox_model_sig <- coxph(
  Surv(time, cid) ~ age + karnof + cd80 + cd820 + symptom + z30 + strata(trt, offtrt),
  data = aids_data
)


ui <- fluidPage(
  titlePanel("AIDS Survival Risk Prediction (Significant Predictors Only)"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("age", "Age (years):", value = 35, min = 0, max = 100),
      numericInput("karnof", "Karnofsky Score (0–100):", value = 90, min = 0, max = 100),
      numericInput("cd80", "CD8 Count at Baseline:", value = 900),
      numericInput("cd820", "CD8 Count at ~20 Weeks:", value = 850),
      selectInput("symptom", "Symptomatic?", choices = c("Asymptomatic" = 0, "Symptomatic" = 1)),
      selectInput("z30", "ZDV in Last 30 Days?", choices = c("No" = 0, "Yes" = 1)),
      selectInput("trt", "Treatment Group:", choices = c("ZDV only" = 0, "ZDV + ddI" = 1, "ZDV + Zal" = 2, "ddI only" = 3)),
      selectInput("offtrt", "Stopped Treatment Early?", choices = c("No" = 0, "Yes" = 1))
    ),
    
    mainPanel(
      h3("Predicted Risk Score (Linear Predictor)"),
      verbatimTextOutput("riskScore"),
      h3("Personalized Survival Curve"),
      plotOutput("survPlot")
    )
  )
)


server <- function(input, output) {
  
  # Reactive expression for new patient data
  newdata_input <- reactive({
    data.frame(
      age = as.numeric(input$age),
      karnof = as.numeric(input$karnof),
      cd80 = as.numeric(input$cd80),
      cd820 = as.numeric(input$cd820),
      symptom = as.numeric(input$symptom),
      z30 = as.numeric(input$z30),
      trt = as.numeric(input$trt),
      offtrt = as.numeric(input$offtrt)
    )
  })
  
  output$riskScore <- renderPrint({
    newdata <- newdata_input()
    lp <- predict(cox_model_sig, newdata = newdata, type = "lp")
    
    cat("Linear Predictor (Risk Score):", round(lp, 3), "\n")
    cat("Interpretation: Higher values indicate increased risk of event (shorter survival).")
  })
  
  output$survPlot <- renderPlot({
    newdata <- newdata_input()
    
    # Generate survival curve for the new patient
    surv_fit <- survfit(cox_model_sig, newdata = newdata)
    
    # Plot using base R or ggplot2 (below is base R)
    plot(
      surv_fit, 
      xlab = "Time", 
      ylab = "Survival Probability", 
      main = "Estimated Survival Curve for Patient",
      col = "blue",
      lwd = 2
    )
    grid()
  })
}


shinyApp(ui = ui, server = server)

