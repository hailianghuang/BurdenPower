#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
require(shinyjs)
source("helper.r", local=TRUE)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   
   # Application title
   titlePanel("Power analysis for burden test"),
   # Sidebar with a slider input for number of bins 
   fluidRow(
     column(3, wellPanel(
       shinyjs::useShinyjs(),
       selectInput("test", "Choose your analysis", c("De novo", "Case/Control")),
#       conditionalPanel(
 #        condition = "input.test == 'Case/Control' | input.simulation ",
         numericInput("N_rep", "Number of replications", 5),
#         ),
       numericInput("N", "Sample size", 5000),
       numericInput("r", "Case:Control ratio", 1, min=0.1, max=10, step=1)),
       actionButton("goButton", "Run analysis"),
	     downloadButton("downloadPDF", "Save figure")
     ),
     column(3, 
       uiOutput("ui")
     ),
     column(6, 
        plotOutput("plot")
     )
  ),
  p("This power calculator was designed by Drs. Hailiang Huang, Stephan Sanders and Benjamin Neale.  Hailiang Huang and Benjamin Neale are from the Analytic and Translational Genetics Unit, Massachusetts General Hospital and the Broad Institute.  Stephan Sanders is from the Department of Psychiatry, UCSF School of Medicine. Source code, documents and contact info available at", a(href="https://github.com/hailianghuang/BurdenPower", "https://github.com/hailianghuang/BurdenPower"))
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
	
  observe({
	
    output$ui <- renderUI({
      if (is.null(input$test))
        return()
      switch(input$test,
             "Case/Control" = wellPanel(
                        numericInput("K", "Disease prevalence", 0.01, min=0.00001, max=0.2, step=0.01),
                        numericInput("s", "Number of variants in gene", 50, min=1, max=10000, step=1),                
                        numericInput("AF_bar", "Average minor allele frequency", 0.001, min=0.00001, max=0.05, step=0.001), 
                        numericInput("R", "Relative risk for functional alleles", 1.5, min=1, max=20, step=0.5), 
                        numericInput("P_cut_single", "p-value threshold for single variant test", 1E-10, min=1E-12, max=0.05, step=1E-10),
                        numericInput("p_cut_burden", "p-value threshold for burden test", 1E-10, min=1E-12, max=0.05, step=1E-10)
             ),
             "De novo" = wellPanel(
                        numericInput("denovo_genome", "Number of de novo mutations per genome", 60, min=10, max=200, step=5), 
                        numericInput("denovo_select", "Number of de novo mutations in selected regions", 20, min=1, max=200, step=5), 
                        numericInput("f_func_genome", "Proportion of functional alleles in genome ", 0.01, min=0.001, max=1, step=0.01), 
                        numericInput("f_func_select", "Proportion of functional alleles in selected regions", 0.03, min=0.001, max=1, step=0.01),
                        numericInput("p_cut_burden.denovo", "p-value threshold", 1E-10, min=1E-12, max=0.05, step=1E-10),
                        checkboxInput("simulation", "add simulation", value=FALSE)
                        
             )
      )
      
    })
    
    observeEvent(input$simulation, {
      if( input$test=="De novo" & !input$simulation ){
        shinyjs::disable("N_rep")
      }else{
        shinyjs::enable("N_rep")
      }      
    })
    
    observeEvent(input$test, {
      if( input$test=="De novo" ){
        shinyjs::disable("N_rep")
        updateNumericInput(session, "N", value=5000)
      }else{
        shinyjs::enable("N_rep")
        updateNumericInput(session, "N", value=100000)
      }      
    })
    
    set.seed(1)
  
    plotRegular <- function(){
      
      n_function <- c(1:input$s)
      
      AF_raw <- rexp(input$s*input$N_rep, 1/input$AF_bar)
      AF_raw <- matrix(AF_raw, input$s, input$N_rep)
      
      AF_save <- lapply(c(1:input$N_rep), function(i_rep, AF_raw){ 
        prob <- 1/AF_raw[,i_rep];
        prob <- prob/sum(prob);
        ret <- sample(AF_raw[,i_rep], prob=prob)
        ret}, AF_raw)
      
      AF_save <- do.call("cbind", AF_save)
      
      power_burden <- c()
      power_single <- c()
      sum_var <- apply(AF_save, 2, function(x){ sum(x*(1-x) ) })
      
      withProgress(message = 'Making plot', value=0, {
        for (i_rep in c(1:input$N_rep)){
          power_burden <- cbind(power_burden, unlist(lapply(n_function, getBurdenPower, AF_save[,i_rep], sum_var[i_rep], input$R, input$K, input$N, input$r, input$p_cut_burden )))
          power_single <- cbind(power_single, unlist(lapply(n_function, getSingleVarPower, AF_save[,i_rep], input$R, input$K, input$N, input$r, input$P_cut_single, i_rep)))
        }
      })
      
      plot(0, 0, type="n", xlab="# of function variants",  ylab="Power", xlim=range(n_function), ylim=c(0,1))
      legend("topleft", c("Burden","Single Variant"), bty="n", lty=c("solid","dashed"), lwd=3, col=c("blue","red"))
      for (i_rep in c(1:input$N_rep)){
        lines(n_function, power_burden[, i_rep], lwd=3,col="blue")
        lines(n_function, power_single[, i_rep], lwd=3, col="red", lty="dashed")
      }
      
    }

    plotDenovo <- function(){
      
      plot(0, 0, type="n", xlab="Relative Risk", ylab="Power", xlim=c(1,10), ylim=c(0,1))
      rr <- seq(1, 10, by=.5)
      power_para_genome <- getDenovoPower_parametric(rr, input$denovo_genome, input$f_func_genome,  input$N, input$r, input$p_cut_burden.denovo)
      power_para_selected <- getDenovoPower_parametric(rr, input$denovo_select , input$f_func_select, input$N, input$r,  input$p_cut_burden.denovo)
      
      lines(rr, power_para_genome, lwd=3,col="blue")
      lines(rr, power_para_selected, lwd=3, col="red")
      
      if(input$simulation){
        withProgress(message = 'Making plot', value=0, {
       
          power_genome <- c()
          power_selected <- c()
          for (i_rep in c(1:input$N_rep)){
            power_genome <- cbind(power_genome,  getDenovoPower(rr, input$denovo_genome, input$f_func_genome,  input$N, input$r, input$p_cut_burden.denovo, i_rep) )
            power_selected <- cbind(power_selected,  getDenovoPower(rr, input$denovo_select , input$f_func_select, input$N, input$r,  input$p_cut_burden.denovo, i_rep) )
          }
        })
      
        for (i_rep in c(1:input$N_rep)){
          lines(rr, power_genome[, i_rep], lwd=1,col="blue", lty="dashed")
          lines(rr, power_selected[, i_rep], lwd=1, col="red", lty="dashed")
        }
        legend("bottomright",c("Full genome (parametric)","Selected regions (parametric)", "Full genome (simulation)","Selected regions (simulation)"),bty="n", lty=c("solid","solid", "dashed", "dashed"),lwd=3,col=c("blue","red","blue","red"))
      }else{
        legend("bottomright",c("Full genome (parametric)","Selected regions (parametric)"),bty="n", lty=c("solid","solid"),lwd=3,col=c("blue","red"))
        
      }
    }
    
    myFigure <- eventReactive(input$goButton, {
      switch(input$test,
             "De novo" =   plotDenovo(),
             "Case/Control" =   plotRegular()  
      )
    })
    
    output$plot <- renderPlot({
     
      myFigure()
    
    }, height=500, width=500)
 
    output$downloadPDF <- downloadHandler(
      filename = function() {
        paste('power-',  Sys.Date(),  '.pdf', sep='')
      },
      content = function(file){
        pdf(file, width=6, height=6)
        switch(input$test,
               "De novo" =   plotDenovo(),
               "Case/Control" =   plotRegular()  
        )
        dev.off()
      },
      contentType = "application/pdf"
    ) 
	})
})

# Run the application 
shinyApp(ui = ui, server = server)

