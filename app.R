library(shiny)
library(zoo)




sbPanel <- sidebarPanel(
  selectInput("cli.ion", "Ion", c("+H", "-H", "+e", "-e", "+Na"), "+H"),
  textInput("cli.ppm","MS1 ppm", 5),
  textInput("cli.acc","MS2 accept ppm", 3),
  textInput("cli.rej","MS2 reject ppm", 10),
  textInput("cli.ff","Fuzzy formula", "C0-30H0-60N0-10O0-6Cl0-5"),
  textInput("cli.args","arguments", "oei exist"),
  actionButton("cli.build", "Make CLI")
  )


mPanel <-   tabsetPanel(
  tabPanel(
    "input",
    textInput("cli", "CLI parameters", ""),
    textInput("m","Parent ion",0),
    actionButton("m.find", "Find parent ion"),
    p("MS1 (2 or 3 row with tab or space separator"),
    checkboxInput("ms1.header", "with header row", TRUE),
    tags$textarea(id="ms1", rows=10, cols=50),
    p("MS2 (2 or 3 row with tab or space separator"), 
    checkboxInput("ms2.header", "with header row", TRUE),
    tags$textarea(id="ms2", rows=10, cols=50),
    actionButton("process", "Process")
  ),
  tabPanel(
    "data",
    dataTableOutput("ms1"),
    dataTableOutput("ms2")
  ),
  tabPanel(
    "output",
    textOutput("cli"),
    dataTableOutput("formulae"),
    selectInput("formula", "Formula",c()),
    dataTableOutput("msms")
  ),
  tabPanel(
    "msms",
    plotOutput("msmsPlot"))
)

ui <- fluidPage(
  title="gsgui 0.0.0.0.0.0.01 alpha",
  sidebarLayout(
    sbPanel,mainPanel(mPanel))

  )

server <- function(input, output, session)
  {  
  r.ms1 <- reactive({
    r.ms1.plain()
    input$m

    ms1.t <- r.ms1.plain()
    m <- as.numeric(input$m)
    ms1.t <- ms1.t[(as.numeric(ms1.t[,1]) > (m -1)) & (as.numeric(ms1.t[,1]) < (m+ 10)),c(1,2)]
  })
  
  r.ms1.plain <- reactive({
    input$ms1    
    input$ms1.header
    
    ms1 <- input$ms1
    ms1.t <- read.table(text=ms1, sep="", header=input$ms1.header)
    ms1.t
  })
  
  r.ms2 <- reactive({
    input$ms2
    input$ms2.header
    
    ms2 <- input$ms2
    
    
    ms2.t <- read.table(text=ms2, sep="", header=input$ms2.header)
    ms2.t <- ms2.t[,c(1,2)]
    ms2.t
  })
  
  pool <- reactiveValues()
  pool$cli <- ""
  
  gfresult <- reactive({
    input$process
    
    if(input$process > 0)
    {
      isolate({
        #try({
          # format the output values
          progress <- shiny::Progress$new(session, min=1, max=15)
          on.exit(progress$close())
          progress$set(message="Processing...")
          cli <- input$cli
          m <- input$m
          ms1.t <- r.ms1()
          ms1.temp <- tempfile("ms1-")
          write.table(ms1.t, file=ms1.temp, sep="\t", row.names=FALSE, col.names=FALSE)
          
          ms2.t <- r.ms2()
          ms2.temp <- tempfile("ms2-")
          write.table(ms2.t, file=ms2.temp, sep="\t", row.names=FALSE, col.names=FALSE)
          # call GenForm
          cl <- sprintf(
            "%s ms=%s msms=%s %s out analyze",
            gfpath, ms1.temp, ms2.temp, cli)
          r1 <- system(cl, intern=TRUE)
          cl <- sprintf(
            "%s ms=%s msms=%s %s out analyze loss",
            gfpath, ms1.temp, ms2.temp, cli)
          r2 <- system(cl, intern=TRUE)
          
          pool$cli <- cl
          return(list(r1=r1, r2=r2))
        #})
      })
    }
    else(return(NA))
  })
  
  output$ms1 <- renderDataTable(r.ms1())
  output$ms2 <- renderDataTable(r.ms2())
  
  gfprocessed <- reactive({
    d <- gfresult()
    
    if(is.na(d))
      return(NA)
    #save(d, file="d.RData")
    d.processed <- lapply(d, function(r)
    {
      # cut the first 4 lines and the last line out
      r<-r[5:(length(r)-1)]
      # all lines not starting with whitespace are a new formula
      lines.nf <- which( substring(r,1,1) != " ")
      # assign all other lines to a formula
      assigned <- unlist(lapply(1:length(r), function(n) max(c(0,lines.nf)[c(0,lines.nf) < n], na.rm=FALSE)))
      formulae <- r[lines.nf]
      msms <- r[-lines.nf]
      msms.split <- split(msms, assigned[-lines.nf])
      formulae.t <- read.table(text=formulae, header=FALSE, sep="\t")
      colnames(formulae.t) <- c("formula", "dppm", "msmv", "msmsmv", "cmv")
      rownames(formulae.t) <- formulae.t$formula
      formulae.t$lineindex <- as.character(lines.nf)
      msms.split.t <- lapply(msms.split, function(msms.formula)
        {
        msms.f.t <- read.table(text=msms.formula, header=FALSE, sep="\t")
        colnames(msms.f.t) <- c("mz", "formula", "dppm")
        msms.f.t[,"mz"] <- na.locf(msms.f.t[,"mz"])
        msms.f.t
      })
      names(msms.split.t) <- names(msms.split)
      return(list(formulae = formulae.t, msms = msms.split.t))
    })
    d.formulae <- d.processed[[1]]$formulae
    d.formulae$sign <- ifelse(sign(d.formulae$dppm)==1, "+", "-")
    d.formulae$dppm <- abs(d.formulae$dppm)
    d.formulae <- d.formulae[,c("formula", "sign", "dppm", "msmv", "msmsmv", "cmv", "lineindex")]
    d.msms <- lapply(1:length(d.processed[[1]]$msms), function(n)
    {
      r <- cbind(d.processed[[1]]$msms[[n]], d.processed[[2]]$msms[[n]])
      if(!all(r[,1] == r[,4])) stop("Error in merging MSMS tables")
      if(!all(r[,3] == r[,6])) stop("Error in merging MSMS tables")
      r <- r[,c(1,2,5,3)]
      colnames(r) <- c("mz", "formula", "loss", "dppm")
      r[,"sign"] <- ifelse(sign(r$dppm)==1, "+", "-")
      r$dppm <- abs(r$dppm)
      r <- r[,c("mz", "formula", "loss", "sign", "dppm")]
      r
    })
    names(d.msms) <- names(d.processed[[1]]$msms)
    list(formulae = d.formulae, msms=d.msms)
  })
  
  
  observe({
    input$m.find
    
    isolate({
      ms1 <- r.ms1.plain()
      m <- ms1[which.max(ms1[,2]),1]
      updateTextInput(session,"m", value=m)
    })
  })
  
  output$formulae <- renderDataTable(gfprocessed()$formulae[,c("formula", "sign", "dppm", "msmv", "msmsmv", "cmv")])
  
  observe({
    d <- gfprocessed()
    if(!is.na(d))
    {
      l <- names(d$msms)
      names(l) <- d$formulae[match(l,d$formulae$lineindex),"formula"]
      updateSelectInput(session,"formula",choices=l, selected=l[[1]])
    }
  })
  
  cli <- reactive({
    input[["cli.build"]]
    isolate({
      #"ion=+H ff=C0-30H0-60N0-10O0-6Cl0-5 ppm=5 acc=3 rej=10"
      sprintf("ion=%s ff=%s ppm=%s acc=%s rej=%s %s",
              input[["cli.ion"]],
              input[["cli.ff"]],
              input[["cli.ppm"]],
              input[["cli.acc"]],
              input[["cli.rej"]],
              input[["cli.args"]])
    })
  })
  
  observe({
    d <- cli()
    updateTextInput(session, "cli", value=d)
  })
  
  msmsTable <- reactive({
    input$formula
    
    isolate({
      d <- gfprocessed()
      d$msms[[input$formula]]
    })
  })
  
  output$msms <- renderDataTable(msmsTable())
  
  output$cli <- renderText(pool$cli)
  #output$gfresults <- renderText(gfnice())
  
  output$msmsPlot <- renderPlot({
    ms2 <- r.ms2()
    msms <- msmsTable()
    
    plot(ms2[,2] ~ ms2[,1], data=ms2, type='h', col="black")
    abline(v=msms$mz, col="green")
    lines(ms2[,2] ~ ms2[,1], data=ms2, type='h', col="black")
  })


}

if(!exists("gfpath")){
  message("-----------------------------------------------------")
  message("Path to GenForm.exe is not set!")
  message("Before running this shiny app, set the GenForm path like this:")
  message("gfpath <- \"C:/Software/GenForm/GenForm.exe\" ")
  #gfpath <- "C:/Software/GenForm/GenForm.exe"
}
message("-----------------------------------------------------")
message("Ignore all warnings and errors here. They are normal.")
message("-----------------------------------------------------")
shinyApp(server=server, ui=ui)


