# Load data ----
library(ggplot2)
library(matrixStats)
library(readr)
header   <- read_delim("GENENAME_Mean.csv", n_max = 1, delim = ",")
muscle <- read_delim('GENENAME_Mean.csv', delim = ",",
                       col_types=paste(c('c', rep('n', (ncol(header)-1))), collapse=""))
muscle <- data.frame(muscle, row.names = 1)
res <- muscle[1:(ncol(muscle)-7)]


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
  
  sidebarLayout(
      sidebarPanel(width = 4,
      h4(textInput("gene1", "Official Gene Name:", value = "MYH2")),
      h4(helpText("Expression relative to median (log2)")),
      tableOutput("means"),
      downloadButton("downloadPlot", "Save plot"),
      downloadButton("downloadData", "Save data")
      ),
      
  mainPanel(
    plotOutput("geneName")
    )
  )
)


# Define server logic ----
server <- function(input, output) {

  select <- reactive({genename <- toupper(input$gene1)
  x      <- c(rep('A1', length(grep('Biopsy',     colnames(muscle)))),  #list of sample types
              rep('A2', length(grep('Blood',      colnames(muscle)))),
              rep('A3', length(grep('Endo',       colnames(muscle)))),
              rep('A4', length(grep('Fiber',      colnames(muscle)))),
              rep('A5', length(grep('Fibroblast', colnames(muscle)))),
              rep('A6', length(grep('Macrophage',        colnames(muscle)))),
              rep('A7', length(grep('Myocyte',    colnames(muscle))))) 
  data   <- data.frame(x=factor(), y=numeric(), Gene=character(), stringsAsFactors=FALSE) #empty dataframe to collect data
  for( i in 1:length(genename)) { 
    y     <- as.numeric(muscle[genename[i],])            #collect data for gene name i
    datay <- cbind.data.frame(x, y, rep(genename[i])) #create table with x="sample type", y="data", "gene name"
    colnames(datay) <- c("x","y","Gene")              #rename column names to make it possible to rbind later
    data  <- rbind.data.frame(data, datay)            #bind the new gene data at the bottom of the previous one
  }
  data$x <- as.factor(data$x)                         #for a box plot, x should be a factor
  data
  })
  
  stats <- reactive({genename <- toupper(input$gene1)
  mean <- cbind(
    rowMeans(res[genename, grep('Biopsy', colnames(res))], na.rm=T),
    rowMeans(res[genename, grep('Blood', colnames(res))], na.rm=T),
    rowMeans(res[genename, grep('Endo', colnames(res))], na.rm=T),
    rowMeans(res[genename, grep('Fiber', colnames(res))], na.rm=T),
    rowMeans(res[genename, grep('Fibroblast', colnames(res))], na.rm=T),
    rowMeans(res[genename, grep('Macrophage', colnames(res))], na.rm=T),
    rowMeans(res[genename, grep('Myocyte', colnames(res))], na.rm=T))
  Sd <- cbind(
    rowSds(as.matrix(res[genename, grep('Biopsy', colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grep('Blood', colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grep('Endo', colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grep('Fiber', colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grep('Fibroblast', colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grep('Macrophage', colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grep('Myocyte', colnames(res))]), na.rm=T))
  nsize <- cbind(
    rowSums(!is.na(res[genename, grep('Biopsy', colnames(res))])),
    rowSums(!is.na(res[genename, grep('Blood', colnames(res))])),
    rowSums(!is.na(res[genename, grep('Endo', colnames(res))])),
    rowSums(!is.na(res[genename, grep('Fiber', colnames(res))])),
    rowSums(!is.na(res[genename, grep('Fibroblast', colnames(res))])),
    rowSums(!is.na(res[genename, grep('Macrophage', colnames(res))])),
    rowSums(!is.na(res[genename, grep('Myocyte', colnames(res))])))
  stats <- data.frame(t(mean), t(Sd), t(nsize))
  stats[,3] <- as.factor(stats[,3])
  colnames(stats) <- c('Mean', 'Sd', 'n')
  rownames(stats) <- c("Muscle\nBiopsy", "Whole Blood", "Primary\nEndothelium",
                       "Isolated\nMuscle Fiber", "Fibroblast", "Blood-derived\nMacrophage", "Primary\nMyocyte")
  stats
  })
  
  plotInput <- function(){
    data <- select()
    ggplot(data, aes(x=x, y=y, fill=x)) +  geom_boxplot() +
      scale_x_discrete(breaks=c("A1","A2","A3","A4","A5","A6", "A7"),
                       labels=c("Muscle\nBiopsy", "Whole Blood", "Primary\nEndothelium",
                                "Isolated\nMuscle Fiber", "Fibroblast", "Blood-derived\nMacrophage", "Primary\nMyocyte")) +
      labs(x="",
           y=paste(toupper(input$gene1), "expression (log2)"),
           title="",
           caption = "www.nicopillon.com") +
      ylim(min(res, na.rm=T), max(res, na.rm=T)) +
      theme(axis.text.x = element_text(face="bold", color="black", size=14, angle=45, hjust=1),
            axis.text.y = element_text(color="black", size=12, angle=0),
            axis.title  = element_text(face="bold", color="black", size=14, angle=0),
            legend.position="none", legend.title = element_blank()) +
      geom_hline(aes(yintercept=0), linetype="dashed", show.legend=F, color="gray50") +
      geom_text(aes(x=8, y=0.35),label="median", hjust=1, size=4, color="gray40") +
      scale_color_manual(values=c("#D55E00", "#CC79A7", "#0072B2", "#D55E00", "#009E73","#E69F00","#D55E00")) +
      scale_fill_manual(values=c("#D55E00", "#CC79A7", "#0072B2", "#D55E00", "#009E73","#E69F00","#D55E00"))
    }
  
  output$geneName <- renderPlot({
    plotInput()
  })
  
  output$means <- renderTable(align="c", spacing="xs", rownames=T, colnames=T,{
    stats()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$gene1, '.jpeg', sep='') },
    content = function(file) {
      png(file,
          units="cm", width=20, height=12, 
          pointsize=12, res=300)
      print(plotInput())
      dev.off()
    })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$gene1, ".csv", sep="") },
    content = function(file) {
      write.csv(stats(), file)
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)