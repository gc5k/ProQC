#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")
#library(BiocManager)
if(length(grep("apple",sessionInfo()$platform, ignore.case = TRUE))>0) {
  system("git rev-list head --max-count 1 > gitTag.txt")
}

gTag=read.table("gitTag.txt")

library(tsne)
library(shiny)
library(preprocessCore)
options(shiny.maxRequestSize=200*1024^2, shiny.launch.browser=T)
source("qcFun.R")
pcXlim=8
pcYlim=8
# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("SWATH proQC [Santa]"),
   h6(paste("Git commit:", gTag[1,1])),
   h6(paste("Contact: chenguobo@gmail.com")),
   hr(),

   # Sidebar with a slider input for number of bins
   fluidRow(
     column(3,
        fileInput("proFile",
                  "Protein matrix file", accept = c(".csv", ".txt"))
     ),
     column(3,
        fileInput("indFile",
                  "Sample information file", accept = c(".csv", ".txt"))
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
                    choices=c("as.is", "QN"))
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
          tabPanel("PCA",
            plotOutput("generic")
          ),
          tabPanel("eSpace",
            sidebarPanel(
              sliderInput("x", "xLab", 1, pcXlim, value=c(1,1), step = 1, dragRange = T),
              sliderInput("y", "yLab", 1, pcYlim, value=c(2,2), step = 1, dragRange = T),
              selectInput("label", "Group",
                          choices=c("Tissue", "Year", "MS_ID")
              )
            ),
            mainPanel(
              plotOutput("distPlot"),
              plotOutput("distPlot2")
            )
          ),
          tabPanel("tSNE",
            sidebarPanel(
              selectInput("tSNElabel", "Group",
                          choices=c("Tissue", "Year", "MS_ID")
              ),
              actionButton("runtSNE",
                "tSNE, Go!", icon = icon("refresh"))
              ),
            mainPanel(
              plotOutput("tSNEplot")
            )
          )
        )
      )
   )
)

MakeProMat <- function(pFile, sFile, mCut) {

  print(paste0("Reading protein matrix [", Sys.time()[1], "]"))
  if (substr(pFile, nchar(pFile)-2, nchar(pFile)) == "csv") {
    prot=read.csv(pFile, as.is = T, header = T)
  } else {
    prot=read.table(pFile, as.is = T, header = T)
  }
  print(paste0(ncol(prot)-1, " samples and ", nrow(prot), " proteins."))

  #remove missing
  ms=array(0, dim=ncol(prot)-1)
  for(i in 1:length(ms)) {
    ms[i]=length(which(is.na(prot[,i+1])))
  }
  idxMS=which(ms> (mCut*nrow(prot)))

  if (length(idxMS)> 0) {
    prot=prot[,-(idxMS+1)]
  } else {
    prot=prot[,-1]
  }
  print(paste0("Removed ", length(idxMS), " samples [missing rate > ", mCut, "]"))

  colnames(prot)=tolower(colnames(prot))

  print(paste0("Reading sample information [", Sys.time()[1], "]"))
  if (substr(sFile, nchar(sFile)-2, nchar(sFile)) == "csv") {
    phe0=read.csv(sFile, as.is = T, header = T)
  } else {
    phe0=read.table(sFile, as.is = T, header = T)
  }
  phe1=phe0
  phe1$ID=tolower(phe1$ID)

  comP=intersect(colnames(prot), phe1$ID)

  if (length(comP) > 0) {
    prot=prot[ ,colnames(prot) %in% comP]
    phe1=phe1[phe1$ID %in% comP, ]
  } else {
    showNotification("No samples in common.", type = c("error"))
  }

  prot=subset(prot[, order(comP)])
  phe1=subset(phe1[order(comP), ])

  return (list(datP=t(prot), proM=prot, pheM=phe1))
}

proQC <- function(qcM, datM) {
  print(paste("Conducting QC", Sys.time()[1]))

  print(paste0("QC::", qcM))
  if(qcM == "QN") {
    qdat=t(normalize.quantiles(t(datM$datP)))
  } else {
    qdat=datM$datP
  }

  return(qdat)
}

proPCA <- function(exDat) {
  print(paste("Conducting PCA", Sys.time()[1]))
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
    datM = currentProMat()
    if (input$qcMethod == "as.is") {
      return(datM$datP)
    } else {
      return(proQC(input$qcMethod, datM))
    }
  })

  currentPC <- reactive({
    exDat = currentQC()
    proPCA(exDat)
  })

  currentDistMat <- reactive({
    datP = currentQC()
    distP = dist(datP, diag=TRUE)
    return (distP)
  })

  currentTsneMat <- reactive({
    distP = currentDistMat()
    tss = tsne(distP, perplexity = 50)
    return (tss)
  })

  observeEvent(input$runtSNE, {
    output$tSNEplot <- renderPlot({
      withProgress(message="SWATH tSNE:", value=0, {
        n=4
        incProgress(1/n, detail = paste0(" Updating protein matrix ..."))
#        datM=currentProMat()

        incProgress(1/n, detail = paste0(" Calculating Euclidean distance ..."))
#        distP=dist(datM$datP, diag=TRUE)

        incProgress(1/n, detail = paste0(" Conducting tSNE. Very slow ..."))
#        tss=tsne(distP, perplexity = 50)

        datM=currentProMat()
        idxLab=which(tolower(colnames(datM$pheM)) == tolower(input$tSNElabel))
        if (length(idxLab) == 0) {
          showNotification(paste0("No tag '", input$tSNElabel, "' was found"), type = c("error"))
          return
        }

        tColors=rainbow(length(unique(datM$pheM[,idxLab])))
        names(tColors) = unique(datM$pheM[,idxLab])
        tss=currentTsneMat()
        layout(matrix(c(1,1,2), 1, 3))
        plot(tss[,1], tss[,2], xlab="tSNE 1", ylab="tSNE 2", col=tColors[datM$pheM[,idxLab]], bty='n', pch=16, cex=0.5)
        plot(x=NULL, y=NULL, xlim=c(0, 1), ylim=c(0, 1), axes = F, xlab="", ylab="")
        legend("topleft", legend = unique(datM$pheM[,idxLab]), col=tColors[datM$pheM[,idxLab]], pch=16, cex=0.5, bty ='n')
      })
    })
  })

  observeEvent(input$proRun, {
    output$summary <- renderPrint({
      dat=currentProMat();

      print(paste0("QC is '", input$qcMethod, "'"))
      print(paste(nrow(dat$datP), "samples and", ncol(dat$datP), "proteins are remained for analysis."))
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
        if (length(idxLab) == 0) {
          showNotification(paste0("No tag ", input$label), type = c("error"))
          return
        }
      #######Eigenvalue plot
        incProgress(1/n, detail = paste0(" Visualization ..."))

        if (length(names(table(PPP1[,idxLab])))>12) {
          showNotification(paste0("Too many levels: ", length(names(table(PPP1[,idxLab])))), "levels.", type = c("error"))
          return
        }
      ######################pcA for tissue
        layout(matrix(1:(2*ceiling(length(names(table(PPP1[,idxLab])))/2)),
                      2, ceiling(length(names(table(PPP1[,idxLab])))/2), byrow = F))
        xL=input$x[1]
        yL=input$y[1]
        colors=rainbow(length(unique(PPP1[,idxLab])))
        for(i in 1:length(names(table(PPP1[,idxLab])))) {
          spIdx=which(PPP1[,idxLab] == names(table(PPP1[,idxLab]))[i])
          plot(eve[spIdx,xL], eve[spIdx,yL],
               xlab=paste0("PC ", xL), ylab=paste0("PC ", yL), bty='n', axes = T,
               xlim=range(eve[, xL], na.rm = TRUE)*1.1,
               ylim=1.1*range(eve[,yL], na.rm = TRUE),
               col=i, pch=16, cex=0.5)
          points(mean(eve[spIdx, xL]), mean(eve[spIdx, yL]),
                 pch=1, cex=1, lwd=3, col="grey")
          legend("bottomleft", legend = names(table(PPP1[,idxLab]))[i], bty = 'n')
          abline(h=0, v=0, lty=2, col=c("violet", "cyan"))
        }
      })
    })

    output$distPlot2 <- renderPlot(height = 200, {
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

        idxLab=which(tolower(colnames(PPP1)) == tolower(input$label))
        #######Eigenvalue plot
        incProgress(1/n, detail = paste0(" Visualization ..."))

        ######################pcA for tissue
        layout(matrix(1:2, 1, 2))
        xL=input$x[1]
        yL=input$y[1]
        varMat=matrix(0, pcXlim, 2)
        varMat[,1] = eva[1:pcXlim]/sum(eva)
        for(i in 1:pcXlim) {
          an=aov(eve[, i]~as.factor(PPP1[,idxLab]))
          varMat[i,2]=summary(an)[[1]][1,2]/sum(summary(an)[[1]][,2])
        }
        barCol=rep("grey", pcXlim)
        barCol[c(xL, yL)] = c("cyan", "violet")
        barplot(ylab="Raw var", varMat[,1], border = F, beside=T, col=barCol, xlab="eSpace")
        barplot(ylab="Weighted var", varMat[,2], border = F, beside=T, col=barCol, xlab="eSpace")
      })
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
