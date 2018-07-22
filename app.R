# Load data ----
options(shiny.sanitize.errors = TRUE)
library(ggplot2)
library(xlsx)

datalist <- readRDS('data/MuscleAtlas_data.Rds')
statslist <- readRDS('data/MuscleAtlas_statslist.Rds')


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
server <- function(input, output, session) {
  
  plotInput <- function(){
    genename <- toupper(input$gene1)
    data <- datalist[[genename]]
    gg <- ggplot(data) +  geom_boxplot(aes(x=x, y=y, fill=x)) +
      scale_x_discrete(breaks=c("A1","A2","A3","A4","A5","A6", "A7", "A8"),
                       labels=c("Whole Blood", "Primary\nEndothelial Cell", "Fibroblast", 
                                "Blood-derived\nMacrophage",
                                "Muscle\nBiopsy", "Isolated\nMuscle Fiber", 
                                "Primary\nMyocyte", "Primary Smooth\nMuscle Cell")) +
      labs(x="",
           y=paste(toupper(input$gene1), "expression (log2)"),
           title="",
           caption = "www.nicopillon.com") +
      ylim(-4, 8) + scale_y_continuous(breaks = round(seq(-4, 8, by=2),1)) +
      theme(axis.text.x = element_text(face="bold", color="black", size=14, angle=45, hjust=1),
            axis.text.y = element_text(color="black", size=12, angle=0),
            axis.title  = element_text(face="bold", color="black", size=14, angle=0),
            legend.position="none", legend.title = element_blank()) +
      geom_hline(aes(yintercept=0), linetype="dashed", show.legend=F, color="gray60") +
      geom_hline(aes(yintercept=8), linetype="dashed", show.legend=F, color="gray60") +
      geom_hline(aes(yintercept=-4), linetype="dashed", show.legend=F, color="gray60") +
      geom_text(aes(x=9, y=0),label="median", hjust=1, vjust=-0.5, size=4, color="gray60") +
      geom_text(aes(x=9, y=8),label="highest", hjust=1, vjust=1.2, size=4, color="gray60") +
      geom_text(aes(x=9, y=-4),label="lowest", hjust=1, vjust=-0.5, size=4, color="gray60")
    return(gg)
  }


  output$geneName <- renderPlot({
    plotInput()
  })
  
  
  output$means <- renderTable(align="c", spacing="xs", rownames=T, colnames=T,{
    genename <- toupper(input$gene1)
    statslist[[genename]]
  })
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$gene1, '_Muscle_Atlas.jpeg', sep='') },
    content = function(file) {
      png(file,
          units="cm", width=20, height=15, 
          pointsize=12, res=300)
      print(plotInput())
      dev.off()
    })
  
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$gene1, "_Muscle_Atlas.xlsx", sep="") },
    content = function(file) {
      genename <- toupper(input$gene1)
      stats <- statslist[[genename]]
      write.xlsx(stats, file)
    })


  }


# Run the app ----
shinyApp(ui = ui, server = server)