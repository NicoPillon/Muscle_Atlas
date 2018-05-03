# Load data ----
library(ggplot2)
library(matrixStats)
library(readr)
library(xlsx)
header  <- read_delim("GENENAME_Mean.csv", n_max = 1, delim = ",")
muscle  <- read_delim('GENENAME_Mean.csv', delim = ",",
                      col_types=paste(c('c', rep('n', (ncol(header)-1))), collapse=""))
muscle <- data.frame(muscle, row.names = 1)
res <- muscle[1:(ncol(muscle)-7)]
x      <- c(rep('A1', length(grep('Blood',         colnames(res)))),  #list of sample types
            rep('A2', length(grep('Endothelial',   colnames(res)))),
            rep('A3', length(grep('Fibroblast',    colnames(res)))),
            rep('A4', length(grep('Macrophage',    colnames(res)))),
            rep('A5', length(grep('Muscle_Biopsy', colnames(res)))),
            rep('A6', length(grep('Muscle_Fiber',  colnames(res)))),
            rep('A7', length(grep('Myocyte',       colnames(res)))),
            rep('A8', length(grep('Smooth_muscle', colnames(res))))) 

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
    plotOutput("geneName", height="550px")
    )
  )
)


# Define server logic ----
server <- function(input, output) {

  start_time <- Sys.time()
  select <- reactive({genename <- toupper(input$gene1)
  data   <- data.frame(x=factor(), y=numeric(), Gene=character(), stringsAsFactors=FALSE) #empty dataframe to collect data
  for( i in 1:length(genename)) { 
    y     <- as.numeric(res[genename[i],])            #collect data for gene name i
    datay <- cbind.data.frame(x, y, rep(genename[i])) #create table with x="sample type", y="data", "gene name"
    colnames(datay) <- c("x","y","Gene")              #rename column names to make it possible to rbind later
    data  <- rbind.data.frame(data, datay)            #bind the new gene data at the bottom of the previous one
  }
  data$x <- as.factor(data$x)                         #for a box plot, x should be a factor
  data
  })
  
  stats <- reactive({genename <- toupper(input$gene1)
  mean <- cbind(
    rowMeans(res[genename, grepl('Blood',        colnames(res))], na.rm=T),
    rowMeans(res[genename, grepl('Endothelial',  colnames(res))], na.rm=T),
    rowMeans(res[genename, grepl('Fibroblast',   colnames(res))], na.rm=T),
    rowMeans(res[genename, grepl('Macrophage',   colnames(res))], na.rm=T),
    rowMeans(res[genename, grepl('Muscle_Biopsy',colnames(res))], na.rm=T),
    rowMeans(res[genename, grepl('Muscle_Fiber', colnames(res))], na.rm=T),
    rowMeans(res[genename, grepl('Myocyte',      colnames(res))], na.rm=T),
    rowMeans(res[genename, grepl('Smooth_muscle',colnames(res))], na.rm=T))
  Sd <- cbind(
    rowSds(as.matrix(res[genename, grepl('Blood',        colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grepl('Endothelial',  colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grepl('Fibroblast',   colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grepl('Macrophage',   colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grepl('Muscle_Biopsy',colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grepl('Muscle_Fiber', colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grepl('Myocyte',      colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[genename, grepl('Smooth_muscle',colnames(res))]), na.rm=T))
  nsize <- cbind(
    rowSums(!is.na(res[genename, grepl('Blood',        colnames(res))])),
    rowSums(!is.na(res[genename, grepl('Endothelial',  colnames(res))])),
    rowSums(!is.na(res[genename, grepl('Fibroblast',   colnames(res))])),
    rowSums(!is.na(res[genename, grepl('Macrophage',   colnames(res))])),
    rowSums(!is.na(res[genename, grepl('Muscle_Biopsy',colnames(res))])),
    rowSums(!is.na(res[genename, grepl('Muscle_Fiber', colnames(res))])),
    rowSums(!is.na(res[genename, grepl('Myocyte',      colnames(res))])),
    rowSums(!is.na(res[genename, grepl('Smooth_muscle',colnames(res))])))
  stats <- data.frame(t(mean), t(Sd), t(nsize))
  stats[,3] <- as.factor(stats[,3])
  colnames(stats) <- c('Mean', 'Sd', 'n')
  rownames(stats) <- c("Whole Blood", "Endothelial Cell", "Fibroblast", "Macrophage", 
                       "Muscle Biopsy", "Muscle Fiber", "Primary Myocyte", "Smooth Muscle Cell")
  stats
  })
  
  plotInput <- function(){
    data <- select()
    ggplot(data, aes(x=x, y=y, fill=x)) +  geom_boxplot() +
      scale_x_discrete(breaks=c("A1","A2","A3","A4","A5","A6", "A7", "A8"),
                       labels=c("Whole Blood", "Primary\nEndothelial Cell", "Fibroblast", 
                                "Blood-derived\nMacrophage",
                                "Muscle\nBiopsy", "Isolated\nMuscle Fiber", 
                                "Primary\nMyocyte", "Primary Smooth\nMuscle Cell")) +
      labs(x="",
           y=paste(toupper(input$gene1), "expression (log2)"),
           title="",
           caption = "www.nicopillon.com") +
      ylim(min(res, na.rm=T), max(res, na.rm=T)) +
      theme(axis.text.x = element_text(face="bold", color="black", size=14, angle=45, hjust=1),
            axis.text.y = element_text(color="black", size=12, angle=0),
            axis.title  = element_text(face="bold", color="black", size=14, angle=0),
            legend.position="none", legend.title = element_blank()) +
      geom_hline(aes(yintercept=0), linetype="dashed", show.legend=F, color="gray60") +
      geom_hline(aes(yintercept=max(res, na.rm=T)), linetype="dashed", show.legend=F, color="gray60") +
      geom_hline(aes(yintercept=min(res, na.rm=T)), linetype="dashed", show.legend=F, color="gray60") +
      geom_text(aes(x=9, y=0),label="median", hjust=1, size=4, color="gray60") +
      geom_text(aes(x=9, y=max(res, na.rm=T)),label="max", hjust=1, size=4, color="gray60") +
      geom_text(aes(x=9, y=min(res, na.rm=T)),label="min", hjust=1, size=4, color="gray60")
    }
  
  output$geneName <- renderPlot({
    plotInput()
  })
  
  output$means <- renderTable(align="c", spacing="xs", rownames=T, colnames=T,{
    stats()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$gene1, '_Muscle_Atlas.jpeg', sep='') },
    content = function(file) {
      png(file,
          units="cm", width=20, height=12, 
          pointsize=12, res=300)
      print(plotInput())
      dev.off()
    })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$gene1, "_Muscle_Atlas.xlsx", sep="") },
    content = function(file) {
      write.xlsx(stats(), file)
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)