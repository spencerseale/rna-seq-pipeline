if(!require(shiny)) install.packages("shiny", repos = "http://cran.us.r-project.org")
if(!require(shinythemes)) install.packages("shinythemes", repos = "http://cran.us.r-project.org")

library(shiny)
library(DT)
library(shinythemes)

setwd('.')
source("./helper.R", echo = FALSE)

#-------------------three treatments-------------------

if (treatment_num == 3) {
# Define UI ----
#fluidpage adjusts data to users display
ui <- fluidPage(theme = shinytheme("cosmo"),
  #NemaMetrix logo
  img(src = "nema_logo.png", height = 80, width = 300),
  
  #title
  titlePanel(sprintf("RNA-Seq Data Report, %s",Sys.Date())),
  
  #nav bar at the top of page listing comparisons
  navbarPage(
    tabPanel(""),
    tabPanel(
      "Group A vs Group B",
      tabsetPanel(
        tabPanel("Heatmap",
                 sidebarLayout(
                   sidebarPanel(
                     helpText(
                       "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                     ),
                     numericInput(
                       "hnum1",
                       "Number of Genes:",
                       min = 2,
                       max = 10000,
                       value = 30
                     )
                   ),
                   mainPanel(plotOutput("h1", height = "700px"))
                 )),
        navbarMenu(
          "Volcano Plot",
          tabPanel("Figure",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                       ),
                       numericInput(
                         "vpval1",
                         "P-value Cutoff:",
                         min = 10e-16,
                         max = 1,
                         value = 10e-6
                       ),
                       sliderInput(
                         "vlfc1",
                         "Log2 Fold Change Cutoff:",
                         min = 0,
                         max = 10,
                         value = 1
                       )
                     ),
                     mainPanel(plotOutput("v1", height = "700px"))
                   )),
          tabPanel("Data",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                       )
                     ),
                     mainPanel(dataTableOutput("vd1"))
                   ))
        ),
        tabPanel("Gene Ontology",
                 sidebarLayout(
                   sidebarPanel(
                     helpText(
                       "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                     ),
                     radioButtons("gorad1", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                   ),
                   mainPanel(dataTableOutput("go1"))
                 ))
      )
    ),
    tabPanel(
      "Group A vs Group C",
      tabsetPanel(
        tabPanel("Heatmap",
                 sidebarLayout(
                   sidebarPanel(
                     helpText(
                       "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                     ),
                     numericInput(
                       "hnum2",
                       "Number of Genes:",
                       min = 2,
                       max = 10000,
                       value = 30
                     )
                   ),
                   mainPanel(plotOutput("h2", height = "700px"))
                 )),
        navbarMenu(
          "Volcano Plot",
          tabPanel("Figure",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                       ),
                       numericInput(
                         "vpval2",
                         "P-value Cutoff:",
                         min = 10e-16,
                         max = 1,
                         value = 10e-6
                       ),
                       sliderInput(
                         "vlfc2",
                         "Log2 Fold Change Cutoff:",
                         min = 0,
                         max = 10,
                         value = 1
                       )
                     ),
                     mainPanel(plotOutput("v2", height = "700px"))
                   )),
          tabPanel("Data",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                       )
                     ),
                     mainPanel(dataTableOutput("vd2"))
                   ))
        ),
        tabPanel("Gene Ontology",
                 sidebarLayout(
                   sidebarPanel(
                     helpText(
                       "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                     ),
                     radioButtons("gorad2", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                   ),
                   mainPanel(dataTableOutput("go2"))
                 ))
      )
    ),
    tabPanel(
      "Group B vs Group C",
      tabsetPanel(
        tabPanel("Heatmap",
                 sidebarLayout(
                   sidebarPanel(
                     helpText(
                       "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                     ),
                     numericInput(
                       "hnum3",
                       "Number of Genes:",
                       min = 2,
                       max = 10000,
                       value = 30
                     )
                   ),
                   mainPanel(plotOutput("h3", height = "700px"))
                 )),
        navbarMenu(
          "Volcano Plot",
          tabPanel("Figure",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                       ),
                       numericInput(
                         "vpval3",
                         "P-value Cutoff:",
                         min = 10e-16,
                         max = 1,
                         value = 10e-6
                       ),
                       sliderInput(
                         "vlfc3",
                         "Log2 Fold Change Cutoff:",
                         min = 0,
                         max = 10,
                         value = 1
                       )
                     ),
                     mainPanel(plotOutput("v3", height = "700px"))
                   )),
          tabPanel("Data",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                       )
                     ),
                     mainPanel(dataTableOutput("vd3"))
                   ))
        ),
        tabPanel("Gene Ontology",
                 sidebarLayout(
                   sidebarPanel(
                     helpText(
                       "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                     ),
                     radioButtons("gorad3", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                   ),
                   mainPanel(dataTableOutput("go3"))
                 ))
      )
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  #heatmaps
  output$h1 <- renderPlot({
    coolmap(logCPM_A_vs_B[o1[1:input$hnum1],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
  })
  
  output$h2 <- renderPlot({
    coolmap(logCPM_A_vs_C[o2[1:input$hnum2],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
  })
  
  output$h3 <- renderPlot({
    coolmap(logCPM_B_vs_C[o3[1:input$hnum3],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
  })
  
  #volcano plots 
  output$v1 <- renderPlot({
    EnhancedVolcano(A_vs_B_lrt$table, FCcutoff = input$vlfc1, pCutoff = input$vpval1, x = "logFC", y = "PValue", lab=rownames(A_vs_B_lrt$table), title="Group A vs Group B", subtitle="+log2fc favors Group A")
  })
  
  output$v2 <- renderPlot({
    EnhancedVolcano(A_vs_C_lrt$table, FCcutoff = input$vlfc2, pCutoff = input$vpval2, x = "logFC", y = "PValue", lab=rownames(A_vs_C_lrt$table), title="Group A vs Group C", subtitle="+log2fc favors Group A")
  })
  
  output$v3 <- renderPlot({
    EnhancedVolcano(B_vs_C_lrt$table, FCcutoff = input$vlfc3, pCutoff = input$vpval3, x = "logFC", y = "PValue", lab=rownames(B_vs_C_lrt$table), title="Group B vs Group C", subtitle="+log2fc favors Group B")
  })
  
  #volcano plot data
  output$vd1 <- renderDataTable({
    datatable(A_vs_B_lrt$table[order(A_vs_B_lrt$table$PValue),])
  })
  
  output$vd2 <- renderDataTable({
    datatable(A_vs_C_lrt$table[order(A_vs_C_lrt$table$PValue),])
  }) 
  
  output$vd3 <- renderDataTable({
    datatable(B_vs_C_lrt$table[order(B_vs_C_lrt$table$PValue),])
  }) 
  
  #gene ontology
  output$go1 <- renderDataTable({
    datatable(topGO(go_A_vs_B, n = 100, sort = input$gorad1))
  })
  
  output$go2 <- renderDataTable({
    datatable(topGO(go_A_vs_C, n = 100, sort = input$gorad2))
  })
  
  output$go3 <- renderDataTable({
    datatable(topGO(go_B_vs_C, n = 100, sort = input$gorad3))
  })
}

}
  
#-------------------two treatments-------------------

if (treatment_num == 2) {
    # Define UI ----
    #fluidpage adjusts data to users display
    ui <- fluidPage(theme = shinytheme("cosmo"),
      #NemaMetrix logo
      img(src = "nema_logo.png", height = 80, width = 300),
      
      #title
      titlePanel(sprintf("RNA-Seq Data Report, %s",Sys.Date())),
      
      #nav bar at the top of page listing comparisons
      navbarPage(
        tabPanel(""),
        tabPanel(
          "Group A vs Group B",
          tabsetPanel(
            tabPanel("Heatmap",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                         ),
                         numericInput(
                           "hnum1",
                           "Number of Genes:",
                           min = 2,
                           max = 10000,
                           value = 30
                         )
                       ),
                       mainPanel(plotOutput("h1", height = "700px"))
                     )),
            navbarMenu(
              "Volcano Plot",
              tabPanel("Figure",
                       sidebarLayout(
                         sidebarPanel(
                           helpText(
                             "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                           ),
                           numericInput(
                             "vpval1",
                             "P-value Cutoff:",
                             min = 10e-16,
                             max = 1,
                             value = 10e-6
                           ),
                           sliderInput(
                             "vlfc1",
                             "Log2 Fold Change Cutoff:",
                             min = 0,
                             max = 10,
                             value = 1
                           )
                         ),
                         mainPanel(plotOutput("v1", height = "700px"))
                       )),
              tabPanel("Data",
                       sidebarLayout(
                         sidebarPanel(
                           helpText(
                             "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                           )
                         ),
                         mainPanel(dataTableOutput("vd1"))
                       ))
            ),
            tabPanel("Gene Ontology",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                         ),
                         radioButtons("gorad1", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                       ),
                       mainPanel(dataTableOutput("go1"))
                     ))
          )
        )
      )
    )

    # Define server logic ----
    server <- function(input, output) {
      
      #heatmaps
      output$h1 <- renderPlot({
        coolmap(logCPM_A_vs_B[o1[1:input$hnum1],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
      })
  
      #volcano plots 
      output$v1 <- renderPlot({
        EnhancedVolcano(A_vs_B_lrt$table, FCcutoff = input$vlfc1, pCutoff = input$vpval1, x = "logFC", y = "PValue", lab=rownames(A_vs_B_lrt$table), title="Group A vs Group B", subtitle="+log2fc favors Group A")
      })
  
      #volcano plot data
      output$vd1 <- renderDataTable({
        datatable(A_vs_B_lrt$table[order(A_vs_B_lrt$table$PValue),])
      })
      
      #gene ontology
      output$go1 <- renderDataTable({
        datatable(topGO(go_A_vs_B, n = 100, sort = input$gorad1))
      })
    }
    
}

#-------------------four treatments-------------------

if (treatment_num == 4) {
  # Define UI ----
  #fluidpage adjusts data to users display
  ui <- fluidPage(theme = shinytheme("cosmo"),
    #NemaMetrix logo
    img(
      src = "nema_logo.png",
      height = 80,
      width = 300
    ),
    
    #title
    titlePanel(sprintf("RNA-Seq Data Report, %s", Sys.Date())),
    
    #nav bar at the top of page listing comparisons
    navbarPage(
      tabPanel(""),
      tabPanel(
        "Group A vs Group B",
        tabsetPanel(
          tabPanel("Heatmap",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                       ),
                       numericInput(
                         "hnum1",
                         "Number of Genes:",
                         min = 2,
                         max = 10000,
                         value = 30
                       )
                     ),
                     mainPanel(plotOutput("h1", height = "700px"))
                   )),
          navbarMenu(
            "Volcano Plot",
            tabPanel("Figure",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                         ),
                         numericInput(
                           "vpval1",
                           "P-value Cutoff:",
                           min = 10e-16,
                           max = 1,
                           value = 10e-6
                         ),
                         sliderInput(
                           "vlfc1",
                           "Log2 Fold Change Cutoff:",
                           min = 0,
                           max = 10,
                           value = 1
                         )
                       ),
                       mainPanel(plotOutput("v1", height = "700px"))
                     )),
            tabPanel("Data",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                         )
                       ),
                       mainPanel(dataTableOutput("vd1"))
                     ))
          ),
          tabPanel("Gene Ontology",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                       ),
                       radioButtons("gorad1", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                     ),
                     mainPanel(dataTableOutput("go1"))
                   ))
        )
      ),
      tabPanel(
        "Group A vs Group C",
        tabsetPanel(
          tabPanel("Heatmap",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                       ),
                       numericInput(
                         "hnum2",
                         "Number of Genes:",
                         min = 2,
                         max = 10000,
                         value = 30
                       )
                     ),
                     mainPanel(plotOutput("h2", height = "700px"))
                   )),
          navbarMenu(
            "Volcano Plot",
            tabPanel("Figure",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                         ),
                         numericInput(
                           "vpval2",
                           "P-value Cutoff:",
                           min = 10e-16,
                           max = 1,
                           value = 10e-6
                         ),
                         sliderInput(
                           "vlfc2",
                           "Log2 Fold Change Cutoff:",
                           min = 0,
                           max = 10,
                           value = 1
                         )
                       ),
                       mainPanel(plotOutput("v2", height = "700px"))
                     )),
            tabPanel("Data",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                         )
                       ),
                       mainPanel(dataTableOutput("vd2"))
                     ))
          ),
          tabPanel("Gene Ontology",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                       ),
                       radioButtons("gorad2", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                     ),
                     mainPanel(dataTableOutput("go2"))
                   ))
        )
      ),
      tabPanel(
        "Group B vs Group C",
        tabsetPanel(
          tabPanel("Heatmap",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                       ),
                       numericInput(
                         "hnum3",
                         "Number of Genes:",
                         min = 2,
                         max = 10000,
                         value = 30
                       )
                     ),
                     mainPanel(plotOutput("h3", height = "700px"))
                   )),
          navbarMenu(
            "Volcano Plot",
            tabPanel("Figure",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                         ),
                         numericInput(
                           "vpval3",
                           "P-value Cutoff:",
                           min = 10e-16,
                           max = 1,
                           value = 10e-6
                         ),
                         sliderInput(
                           "vlfc3",
                           "Log2 Fold Change Cutoff:",
                           min = 0,
                           max = 10,
                           value = 1
                         )
                       ),
                       mainPanel(plotOutput("v3", height = "700px"))
                     )),
            tabPanel("Data",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                         )
                       ),
                       mainPanel(dataTableOutput("vd3"))
                     ))
          ),
          tabPanel("Gene Ontology",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                       ),
                       radioButtons("gorad3", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                     ),
                     mainPanel(dataTableOutput("go3"))
                   ))
        )
      ),
      tabPanel(
        "Group A vs Group D",
        tabsetPanel(
          tabPanel("Heatmap",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                       ),
                       numericInput(
                         "hnum4",
                         "Number of Genes:",
                         min = 2,
                         max = 10000,
                         value = 30
                       )
                     ),
                     mainPanel(plotOutput("h4", height = "700px"))
                   )),
          navbarMenu(
            "Volcano Plot",
            tabPanel("Figure",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                         ),
                         numericInput(
                           "vpval4",
                           "P-value Cutoff:",
                           min = 10e-16,
                           max = 1,
                           value = 10e-6
                         ),
                         sliderInput(
                           "vlfc4",
                           "Log2 Fold Change Cutoff:",
                           min = 0,
                           max = 10,
                           value = 1
                         )
                       ),
                       mainPanel(plotOutput("v4", height = "700px"))
                     )),
            tabPanel("Data",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                         )
                       ),
                       mainPanel(dataTableOutput("vd4"))
                     ))
          ),
          tabPanel("Gene Ontology",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                       ),
                       radioButtons("gorad4", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                     ),
                     mainPanel(dataTableOutput("go4"))
                   ))
        )
      ),
      tabPanel(
        "Group B vs Group D",
        tabsetPanel(
          tabPanel("Heatmap",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                       ),
                       numericInput(
                         "hnum5",
                         "Number of Genes:",
                         min = 2,
                         max = 10000,
                         value = 30
                       )
                     ),
                     mainPanel(plotOutput("h5", height = "700px"))
                   )),
          navbarMenu(
            "Volcano Plot",
            tabPanel("Figure",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                         ),
                         numericInput(
                           "vpval5",
                           "P-value Cutoff:",
                           min = 10e-16,
                           max = 1,
                           value = 10e-6
                         ),
                         sliderInput(
                           "vlfc5",
                           "Log2 Fold Change Cutoff:",
                           min = 0,
                           max = 10,
                           value = 1
                         )
                       ),
                       mainPanel(plotOutput("v5", height = "700px"))
                     )),
            tabPanel("Data",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                         )
                       ),
                       mainPanel(dataTableOutput("vd5"))
                     ))
          ),
          tabPanel("Gene Ontology",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                       ),
                       radioButtons("gorad5", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                     ),
                     mainPanel(dataTableOutput("go5"))
                   ))
        )
      ),
      tabPanel(
        "Group C vs Group D",
        tabsetPanel(
          tabPanel("Heatmap",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "In the text box below, input the number of genes to be included in the heatmap. The genes are added to the figure in the order of most significantly differentially expressed. Genes are grouped via hierarchical clustering. Note: as more genes are added to the figure the gene names stack and are not reliably tied to their associated panel (~50 genes)."
                       ),
                       numericInput(
                         "hnum6",
                         "Number of Genes:",
                         min = 2,
                         max = 10000,
                         value = 30
                       )
                     ),
                     mainPanel(plotOutput("h6", height = "700px"))
                   )),
          navbarMenu(
            "Volcano Plot",
            tabPanel("Figure",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "Adjust the location of the dashed cutoff lines with the controls below. Genes surpassing the cutoffs will automatically be labeled without overcrowding. However, this may cause some genes to not be labeled in the figure and thus the supporting data table can be used to identify these genes."
                         ),
                         numericInput(
                           "vpval6",
                           "P-value Cutoff:",
                           min = 10e-16,
                           max = 1,
                           value = 10e-6
                         ),
                         sliderInput(
                           "vlfc6",
                           "Log2 Fold Change Cutoff:",
                           min = 0,
                           max = 10,
                           value = 1
                         )
                       ),
                       mainPanel(plotOutput("v6", height = "700px"))
                     )),
            tabPanel("Data",
                     sidebarLayout(
                       sidebarPanel(
                         helpText(
                           "This table contains the quantitative data that the volcano plot is built on. All genes listed in the table are ordered by p-value. Use the search bar to search for genes of interest."
                         )
                       ),
                       mainPanel(dataTableOutput("vd6"))
                     ))
          ),
          tabPanel("Gene Ontology",
                   sidebarLayout(
                     sidebarPanel(
                       helpText(
                         "This table provides the number of differentially expressed genes associated with each Gene Ontology term. The numbers in the 'Up' column are the number of genes overexpressed in the first group relative to the other group in the selected pairwise comparison (i.e. Group A vs. Group B, Group A is referred to as 'Up'). Choose below whether the table should be ordered by lowest 'Up' p-value or lowest 'Down' p-value."
                       ),
                       radioButtons("gorad6", "Order GO Terms By:", c("Up" = "up", "Down" = "down"))
                     ),
                     mainPanel(dataTableOutput("go6"))
                   ))
        )
      )
      
    )
  )
  
  # Define server logic ----
  server <- function(input, output) {
    
    #heatmaps
    output$h1 <- renderPlot({
      coolmap(logCPM_A_vs_B[o1[1:input$hnum1],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
    })
    
    output$h2 <- renderPlot({
      coolmap(logCPM_A_vs_C[o2[1:input$hnum2],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
    })
    
    output$h3 <- renderPlot({
      coolmap(logCPM_B_vs_C[o3[1:input$hnum3],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
    })
    
    output$h4 <- renderPlot({
      coolmap(logCPM_A_vs_D[o4[1:input$hnum4],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
    })
    
    output$h5 <- renderPlot({
      coolmap(logCPM_B_vs_D[o5[1:input$hnum5],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
    })
    
    output$h6 <- renderPlot({
      coolmap(logCPM_C_vs_D[o6[1:input$hnum6],], show.dendrogram = "none", margins = c(8,8), cexRow=0.75, cexCol=1, srtCol=45, cluster.by = "expression level", col="redblue")
    })
        
    #volcano plots 
    output$v1 <- renderPlot({
      EnhancedVolcano(A_vs_B_lrt$table, FCcutoff = input$vlfc1, pCutoff = input$vpval1, x = "logFC", y = "PValue", lab=rownames(A_vs_B_lrt$table), title="Group A vs Group B", subtitle="+log2fc favors Group A")
    })
    
    output$v2 <- renderPlot({
      EnhancedVolcano(A_vs_C_lrt$table, FCcutoff = input$vlfc2, pCutoff = input$vpval2, x = "logFC", y = "PValue", lab=rownames(A_vs_C_lrt$table), title="Group A vs Group C", subtitle="+log2fc favors Group A")
    })
    
    output$v3 <- renderPlot({
      EnhancedVolcano(B_vs_C_lrt$table, FCcutoff = input$vlfc3, pCutoff = input$vpval3, x = "logFC", y = "PValue", lab=rownames(B_vs_C_lrt$table), title="Group B vs Group C", subtitle="+log2fc favors Group B")
    })
    
    output$v4 <- renderPlot({
      EnhancedVolcano(A_vs_D_lrt$table, FCcutoff = input$vlfc4, pCutoff = input$vpval4, x = "logFC", y = "PValue", lab=rownames(A_vs_D_lrt$table), title="Group A vs Group D", subtitle="+log2fc favors Group A")
    })
    
    output$v5 <- renderPlot({
      EnhancedVolcano(B_vs_D_lrt$table, FCcutoff = input$vlfc5, pCutoff = input$vpval5, x = "logFC", y = "PValue", lab=rownames(B_vs_D_lrt$table), title="Group B vs Group D", subtitle="+log2fc favors Group B")
    })
    
    output$v6 <- renderPlot({
      EnhancedVolcano(C_vs_D_lrt$table, FCcutoff = input$vlfc6, pCutoff = input$vpval6, x = "logFC", y = "PValue", lab=rownames(C_vs_D_lrt$table), title="Group C vs Group D", subtitle="+log2fc favors Group C")
    })
    
    #volcano plot data
    output$vd1 <- renderDataTable({
      datatable(A_vs_B_lrt$table[order(A_vs_B_lrt$table$PValue),])
    })
    
    output$vd2 <- renderDataTable({
      datatable(A_vs_C_lrt$table[order(A_vs_C_lrt$table$PValue),])
    }) 
    
    output$vd3 <- renderDataTable({
      datatable(B_vs_C_lrt$table[order(B_vs_C_lrt$table$PValue),])
    }) 
    
    output$vd4 <- renderDataTable({
      datatable(A_vs_D_lrt$table[order(A_vs_D_lrt$table$PValue),])
    }) 
    
    output$vd5 <- renderDataTable({
      datatable(B_vs_D_lrt$table[order(B_vs_D_lrt$table$PValue),])
    }) 
    
    output$vd6 <- renderDataTable({
      datatable(C_vs_D_lrt$table[order(C_vs_D_lrt$table$PValue),])
    }) 
    
    #gene ontology
    output$go1 <- renderDataTable({
      datatable(topGO(go_A_vs_B, n = 100, sort = input$gorad1))
    })
    
    output$go2 <- renderDataTable({
      datatable(topGO(go_A_vs_C, n = 100, sort = input$gorad2))
    })
    
    output$go3 <- renderDataTable({
      datatable(topGO(go_B_vs_C, n = 100, sort = input$gorad3))
    })
    
    output$go4 <- renderDataTable({
      datatable(topGO(go_A_vs_D, n = 100, sort = input$gorad4))
    })
    
    output$go5 <- renderDataTable({
      datatable(topGO(go_B_vs_D, n = 100, sort = input$gorad5))
    })
    
    output$go6 <- renderDataTable({
      datatable(topGO(go_C_vs_D, n = 100, sort = input$gorad6))
    })
  }
  
}

#defining the app UI and server
shinyApp(ui = ui, server = server)