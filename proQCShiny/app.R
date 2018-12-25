#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")

library(shiny)
library(preprocessCore)
options(shiny.maxRequestSize=200*1024^2, shiny.launch.browser=T)
source("qcFun.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("SWATH proQC [Santa Claus]"),

   # Sidebar with a slider input for number of bins
   fluidRow(
     column(3,
        fileInput("proFile",
                  "Protein matrix")
     ),
     column(3,
        fileInput("indFile",
                  "Individal info")
     )
   ),
   hr(),
   fluidRow(
     column(3,
       sliderInput("mRate",
                        "Missing rate",
                        0, 1, 0.85)
     ),
     column(3,
       selectInput("qcMethod", "Quality control",
                    choices=c("as.is", "QN", "Combat", "LM"))
     )
   ),
   fluidRow(
     column(4, offset = 6,
        actionButton("proRun",
                     "Go, Santa Go!", icon = icon("refresh"))
     )
   ),
      # Show a plot of the generated distribution
   fluidRow(
      mainPanel(
        tabsetPanel(type = 'tabs',
          tabPanel("Summary",
            verbatimTextOutput("summary")
          ),
          tabPanel("Generic",
            plotOutput("generic", height = "200px")
          ),
          tabPanel("eSpace",
            sidebarPanel(
              sliderInput("x", "xLab", 1, 8, 1, step = 1),
              sliderInput("y", "yLab", 1, 8, 2, step = 1),
              selectInput("label", "Label",
                          choices=c("Tissue", "Year", "MS_ID")
              )
            ),
            mainPanel(
              plotOutput("distPlot")
            )
          )
        )
      )
   )
)

MakeProMat <- function(pFile, sFile, mCut) {

  print(paste("Making Protein matrix", Sys.time()[1]))
  prot=read.table(pFile, as.is = T, header = T)
  #remove missing

  dat=t(prot[,c(2:ncol(prot))])
  ms=array(0, dim=nrow(dat))
  for(i in 1:length(ms)) {
    ms[i]=length(which(is.na(dat[i,])))
  }
  idxMS=which(ms>(mCut*ncol(dat)))

  dat=dat[-idxMS,]
  prot=prot[,-(idxMS+1)]

  PPP0=read.csv(sFile, as.is = T, header = T)
  PPP1=PPP0
  PPP1$PPPA_ID=tolower(PPP1$PPPA_ID)
  PPP1=PPP1[-idxMS,]
  print(colnames(prot)[1:3])
  return (list(datP=dat, proM=prot, pheM=PPP1))
}

proQC <- function(qcM, datM) {
  print(paste("Concudting QC", Sys.time()[1]))
  prot=datM$proM
  PPP1=datM$pheM

  if(qcM == "QN") {
    print("QC::QN")
    qdat=t(normalize.quantiles(t(datM$datP)))
  } else {
    print("Other QC")
    qdat=datM$datP
  }
  colnames(qdat)=prot$protein_group
  rownames(qdat)=rownames(t(prot[,c(2:ncol(prot))]))

  #line up data
  exDat=matrix(0, nrow(qdat), ncol(qdat))
  rID=matrix(0, nrow(qdat), 1)
  for (i in 1:nrow(PPP1)) {
    idx=which(rownames(qdat) == PPP1$PPPA_ID[i])
    rID[i,1]=idx
    exDat[i,] = qdat[idx,]
  }
  rownames(exDat) = PPP1$PPPA_ID[rID[,1]]
  colnames(exDat)=colnames(qdat)
  return(exDat)
}

proPCA <- function(exDat) {
  print(paste("Concudting PCA", Sys.time()[1]))
  p_CorM=proCorMatrix_c(exDat)
  p_Eg=eigen(p_CorM)
  return(list(pCorM=p_CorM, pEg=p_Eg))
}

# Define server logic required to draw a histogram
server <- function(input, output) {

  currentProMat <- reactive({
    MakeProMat(input$proFile$datapath[1], input$indFile$datapath[1], as.numeric(input$mRate))
  })

  currentQC <- reactive({
    datM=currentProMat()
    if (input$qcMethod == "as.is") {
      return(datM$datP)
    } else {
      return(proQC(input$qcMethod, datM))
    }
  })

  currentPC <- reactive({
    exDat=currentQC()
    proPCA(exDat)
  })

  observeEvent(input$proRun, {
    output$summary <- renderPrint({
      print(paste("Protein matrix file:", input$proFile$name))
      print(paste("Protein matrix file:", input$indFile$name))
      print(paste("Removed individuals missing rate >", as.numeric(input$mRate)))
      print(paste("QC is", input$qcMethod))
    })

    output$generic <- renderPlot({

      withProgress(message="SWATH proQC:", value=0, {
        n=3
        incProgress(1/n, detail = paste0(" Reading ..."))
        incProgress(1/n, detail = paste0(" PCA ..."))

        pc=currentPC()
        grm=pc$pCorM
        eve=pc$pEg$vectors
        eva=pc$pEg$values

        incProgress(1/n, detail = paste0(" Visualizing ..."))

        layout(matrix(1:2, 1, 2))
        hist(grm[col(grm)<row(grm)], breaks=50, main = "Protein relationship matrix", xlab="Score")
        cutEva=eva[which(eva>1)]
        barplot(cutEva, main=paste("Top", length(cutEva), "eigenvalues"), xlab="Eigenvalue", col="grey", border = "NA")
        abline(h=1, col="red", lty=2)
      })
    })

    output$distPlot <- renderPlot({
      ##read data
      withProgress(message="SWATH proQC:", value=0, {
        n = 3
        incProgress(1/n, detail = paste0(" Reading ..."))
        incProgress(1/n, detail = paste0(" PCA ..."))

        pc=currentPC()
        grm=pc$pCorM
        eve=pc$pEg$vectors
        eva=pc$pEg$values
        PPP1=currentProMat()$pheM

        if (input$label == "") {
          Lab="Tissue"
        } else {
          Lab=input$label
        }

        idxLab=which(tolower(colnames(PPP1)) == tolower(input$label))
        print(idxLab)
      #######Eigenvalue plot
        incProgress(1/n, detail = paste0(" Visualization ..."))

      ######################pcA for tissue
        layout(matrix(1:length(names(table(PPP1[,idxLab]))),
                      2, ceiling(length(names(table(PPP1[,idxLab])))/2), byrow = F))
        for(i in 1:length(names(table(PPP1[,idxLab])))) {
          spIdx=which(PPP1[,idxLab] == names(table(PPP1[,idxLab]))[i])
          plot(eve[spIdx,input$x], eve[spIdx,input$y],
               xlab=paste0("PC ", input$x), ylab=paste0("PC ", input$y), bty='n', axes = T,
               xlim=range(eve[,input$x], na.rm = TRUE)*1.1,
               ylim=1.1*range(eve[,input$x], na.rm = TRUE),
               col=i, pch=16, cex=0.5)
          points(mean(eve[spIdx,input$x]), mean(eve[spIdx,input$y]),
                 pch=1, cex=1, lwd=3, col="grey")
          legend("bottomleft", legend = names(table(PPP1[,idxLab]))[i], bty = 'n')
        }
      })
    })

  })
}

# Run the application
shinyApp(ui = ui, server = server)

