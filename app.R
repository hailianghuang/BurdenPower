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
       selectInput("xaxis", "Choose the x-axis", c("Relative risk", "Sample size")),
       numericInput("N_rep", "Number of replications", 5),
       numericInput("N", "Sample size", 5000),
       numericInput("r", "Case:Control ratio", 1, min=0.1, max=10, step=1),
       numericInput("R", "Maximum relative risk", 10, min=1.5, max=10)
      ),
       actionButton("goButton", "Run analysis"),
	     downloadButton("downloadPDF", "Save figure"),
       downloadButton("downloadDat", "Save data")
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
                        numericInput("f_gene", "Number of variants that are functional", 25, min=1, max=10000, step=1), 
                        numericInput("AF_bar", "Average minor allele frequency", 0.001, min=0.00001, max=0.05, step=0.001), 
                        numericInput("P_cut_single", "p-value threshold for single variant test", 1E-10, min=1E-12, max=0.05, step=1E-10),
                        numericInput("p_cut_burden", "p-value threshold for burden test", 1E-10, min=1E-12, max=0.05, step=1E-10),
                        checkboxInput("showCI", "show 95% confidence interval", value=FALSE)
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
    
    observeEvent(input$xaxis, {
      if( input$xaxis=="Relative risk" ){
        if( input$test=="De novo" ){
          updateNumericInput(session, "N", value=5000, label="Sample size")
        }else{
          updateNumericInput(session, "N", value=100000, label="Sample size")
        } 
        updateNumericInput(session, "R", value=10, label="Maximum relative risk")
      }else{
        if( input$test=="De novo" ){
          updateNumericInput(session, "N", value=10000, label="Maximum sample size")
        }else{
          updateNumericInput(session, "N", value=200000, label="Maximum sample size")
        } 
        updateNumericInput(session, "R", value=3, label="Relative risk")
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
  
    plotRegular <- function(download=F, file=NULL){
    
      if(input$f_gene > input$s){
        plot(0, 0, type="n", xlab="",  ylab="", xlim=c(0,1), ylim=c(0,1), axes=F)
        text(0,1, adj=0, "Number of functonal variant must be less than total number of variants.", col="red")
        return()
      }
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
      
      rr <- seq(1, input$R, by=0.5)
      nn <- seq(100, input$N, by=1000)
      withProgress(message = 'Making plot', value=0, {
        for (i_rep in c(1:input$N_rep)){
          incProgress(1/input$N_rep, detail = paste("Replication", i_rep))
          if(input$xaxis=="Relative risk"){
            power_burden <- cbind(power_burden, unlist(lapply(rr, function(R, f, AF, sum_var, K, N, r, p_cut_burden){getBurdenPower(f, AF, sum_var, R, K, N, r, p_cut_burden)}, input$f_gene, AF_save[,i_rep], sum_var[i_rep], input$K, input$N, input$r, input$p_cut_burden )))
            power_single <- cbind(power_single, unlist(lapply(rr, function(R, f, AF, K, N, r, P_cut_single) {getSingleVarPower(f, AF, R, K, N, r, P_cut_single)}, input$f_gene, AF_save[,i_rep],  input$K, input$N, input$r, input$P_cut_single)))
          }else{
            power_burden <- cbind(power_burden, unlist(lapply(nn, function(N, f, AF, sum_var, K, R, r, p_cut_burden){getBurdenPower(f, AF, sum_var, R, K, N, r, p_cut_burden)}, input$f_gene, AF_save[,i_rep], sum_var[i_rep], input$K, input$R, input$r, input$p_cut_burden )))
            power_single <- cbind(power_single, unlist(lapply(nn, function(N, f, AF, K, R, r, P_cut_single) {getSingleVarPower(f, AF, R, K, N, r, P_cut_single)}, input$f_gene, AF_save[,i_rep],  input$K, input$R, input$r, input$P_cut_single)))
          }
          
        }
      })
      
      if(input$xaxis=="Relative risk"){
        xaxis <- rr
        xlabel <- "Relative risk"
      }else{
        xaxis <- nn
        xlabel <- "Sample size"
      }
      
      if(download){
        writeLines(paste("#Analysis:", input$test), file)
        writeLines(paste("#x axis:", input$xaxis), file)
        writeLines(paste("#Sample size:", input$N), file)
        writeLines(paste("#Case/control ratio:", input$r), file)
        writeLines(paste("#Relative risk:", input$R), file)
        writeLines(paste("#Prevalence:", input$K), file)
        writeLines(paste("#Number of variants in gene:", input$s), file)
        writeLines(paste("#Number of variants that are functional:", input$f_gene), file)
        writeLines(paste("#Average minor allele frequency:", input$AF_bar), file)
        writeLines(paste("#p-value threshold for single variant test:", input$P_cut_single), file)
        writeLines(paste("#p-value threshold for burden test:", input$P_cut_burden), file)
        
        if(input$xaxis=="Relative risk"){
          dat <-  data.frame(status=rep(c("burden", "single"), each=length(rr)) ,  RR=rep(rr, 2), rep=rbind(power_burden, power_single))
        }else{
          dat <-  data.frame(status=rep(c("burden", "single"), each=length(nn)) , N=rep(nn, 2), rep=rbind(power_burden, power_single))
        }
        write.table(dat, file, quote=F, col.names=T, row.names=F, sep="\t")
        
      }else{
        
        plot(0, 0, type="n", xlab=xlabel,  ylab="Power", xlim=range(xaxis), ylim=c(0,1))
        
        lines(xaxis, apply(power_burden, 1, mean), lwd=3,col="blue")
        lines(xaxis, apply(power_single, 1, mean), lwd=3,col="red")
        
        if(input$showCI){
          sd <-  apply(power_burden, 1, sd)/sqrt(input$N_rep) * 1.96
          lines(xaxis, pmin(apply(power_burden, 1, mean) +sd,1), lwd=1,col="blue", lty="dashed")
          lines(xaxis, pmax(apply(power_burden, 1, mean) -sd,0) , lwd=1,col="blue", lty="dashed")
          sd <-  apply(power_single, 1, sd)/sqrt(input$N_rep) * 1.96
          lines(xaxis, pmin(apply(power_single, 1, mean) +sd, 1), lwd=1,col="red", lty="dashed")
          lines(xaxis, pmax(apply(power_single, 1, mean) -sd, 0), lwd=1,col="red", lty="dashed")
        }
        legend("topleft", c("Burden","Single Variant"), bty="n", lty=c("solid","solid"), lwd=3, col=c("blue","red"))
        
      }
      
    }

    plotDenovo <- function(download=F, file=NULL){
      
      if(input$f_func_genome * input$denovo_genome < input$denovo_select * input$f_func_select ){
        plot(0, 0, type="n", xlab="",  ylab="", xlim=c(0,1), ylim=c(0,1), axes=F)
        text(0,1, adj=0, "Select a smaller proportion of functional alleles in selected regions", col="red")
        return()
      }
      
      rr <- seq(1, input$R, by=.5)
      nn <- seq(100, input$N, by=1000)
      if(input$xaxis=="Relative risk"){
        power_para_genome <- getDenovoPower_parametric(rr, input$denovo_genome, input$f_func_genome,  input$N, input$r, input$p_cut_burden.denovo)
        power_para_selected <- getDenovoPower_parametric(rr, input$denovo_select , input$f_func_select, input$N, input$r,  input$p_cut_burden.denovo)
      }else{
        power_para_genome <- getDenovoPower_parametric(input$R, input$denovo_genome, input$f_func_genome, nn, input$r, input$p_cut_burden.denovo)
        power_para_selected <- getDenovoPower_parametric(input$R, input$denovo_select , input$f_func_select, nn, input$r,  input$p_cut_burden.denovo)
      }
      
      if(input$simulation){
        withProgress(message = 'Making plot', value=0, {
          
          power_genome <- c()
          power_selected <- c()
          for (i_rep in c(1:input$N_rep)){
            incProgress(1/input$N_rep, detail = paste("Replication", i_rep))
            if(input$xaxis=="Relative risk"){
              power_genome <- cbind(power_genome, unlist(lapply(rr, function(R, q,f,N,r, p_cut_denovo){denovoGenPower(R, q,f,N*r/(1+r),N/(1+r), p_cut_denovo) }, input$denovo_genome,input$f_func_genome,input$N, input$r, input$p_cut_burden.denovo )) )
              power_selected <- cbind(power_selected, unlist(lapply(rr, function(R, q,f,N,r, p_cut_denovo){denovoGenPower(R, q,f,N*r/(1+r),N/(1+r), p_cut_denovo) }, input$denovo_select,input$f_func_select,input$N, input$r, input$p_cut_burden.denovo )) )  
            }else{
              power_genome <- cbind(power_genome, unlist(lapply(nn, function(N, q,f,R,r, p_cut_denovo){denovoGenPower(R, q,f,N*r/(1+r),N/(1+r), p_cut_denovo) }, input$denovo_genome, input$f_func_genome, input$R, input$r, input$p_cut_burden.denovo )) )
              power_selected <- cbind(power_selected, unlist(lapply(nn, function(N, q,f,R,r, p_cut_denovo){denovoGenPower(R, q,f,N*r/(1+r),N/(1+r), p_cut_denovo) }, input$denovo_select, input$f_func_select, input$R, input$r, input$p_cut_burden.denovo )) )  
            }
          }
        })
        
      }
      
      if(input$xaxis=="Relative risk"){
        xaxis <- rr
        xlabel <- "Relative risk"
      }else{
        xaxis <- nn
        xlabel <- "Sample size"
        
      }
      
      if(download){
        writeLines(paste("#Analysis:", input$test), file)
        writeLines(paste("#x axis:", input$xaxis), file)
        writeLines(paste("#Sample size:", input$N), file)
        writeLines(paste("#Case/control ratio:", input$r), file)
        writeLines(paste("#Relative risk:", input$R), file)
        writeLines(paste("#Number of de novo mutations per genome:", input$denovo_genome), file)
        writeLines(paste("#Number of de novo mutations in selected regions:", input$denovo_select), file)
        writeLines(paste("#Proportion of functional alleles in genome:", input$f_func_genome), file)
        writeLines(paste("#Proportion of functional alleles in selected regions:", input$f_func_select), file)
        writeLines(paste("#p-value threshold:", input$p_cut_burden.denovo), file)

        if(input$xaxis=="Relative risk"){
          dat <-  data.frame(status=rep(c("whole_genome", "selected_region"), each=length(rr)) ,  RR=rep(rr, 2), power_parametric=c(power_para_genome, power_para_selected))
        }else{
          dat <-  data.frame(status=rep(c("whole_genome", "selected_region"), each=length(nn)) , N=rep(nn, 2), power_parametric=c(power_para_genome, power_para_selected))
        }
        if(input$simulation){
          dat <- cbind(dat, rep=rbind(power_genome, power_selected))
        }
        write.table(dat, file, quote=F, col.names=T, row.names=F, sep="\t")
      }else{
        plot(0, 0, type="n", xlab=xlabel, ylab="Power", xlim=range(xaxis), ylim=c(0,1))
        lines(xaxis, power_para_genome, lwd=3,col="blue")
        lines(xaxis, power_para_selected, lwd=3, col="red")
        if(input$simulation){
          for (i_rep in c(1:input$N_rep)){
            lines(xaxis, power_genome[, i_rep], lwd=1,col="blue", lty="dashed")
            lines(xaxis, power_selected[, i_rep], lwd=1, col="red", lty="dashed")
          }
          legend("bottomright",c("Full genome (parametric)","Selected regions (parametric)", "Full genome (simulation)","Selected regions (simulation)"),bty="n", lty=c("solid","solid", "dashed", "dashed"),lwd=3,col=c("blue","red","blue","red"))
        }else{
          legend("bottomright",c("Full genome (parametric)","Selected regions (parametric)"),bty="n", lty=c("solid","solid"),lwd=3,col=c("blue","red"))
        }
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
    
    output$downloadDat <- downloadHandler(
      filename = function() {
        paste('power-',  Sys.Date(),  '.txt', sep='')
      },
      content = function(filename){
        con <- file(filename,  open="wt")
        switch(input$test,
               "De novo" =   plotDenovo(TRUE, con),
               "Case/Control" =   plotRegular(TRUE, con)  
        )
        close(con)
      },
      contentType = "application/octet-stream"
    ) 
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

