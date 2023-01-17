
library(shiny)
library(rjson)
library(ggVennDiagram)
library(ggplot2)
library(plotly)
library(DT)
path <- getwd()

# ------ Read in the JSON file to get variants list and constellation ------
#variants <- fromJSON(file = "sgRNA_list.json")
variants <- fromJSON(file = file.path(paste(path,"/sgRNA_list.json",sep = "")))

# ------ UI ------
ui <- shinyUI(fluidPage(

  titlePanel("Intersecting ncsgRNAs of SARS-CoV-2"),
  
  sidebarLayout(
    # --- Input ---
    sidebarPanel(width = 2, 
        checkboxGroupInput("variant", label = "Variant:",
                           choiceNames = names(variants),#selected = c("B.1.1.529 (Omicron)", "B.1.617.2 (Delta)"),
                           choiceValues = names(variants)),
        downloadButton("download1", "Download Graph"),
        downloadButton("download2", "Download Table")
    ),
    # --- Output ---
    mainPanel(width = 10,
      fluidRow(
        plotlyOutput("venn", width = "100%"), style = "height:800px"),
      fluidRow(
        column(
          dataTableOutput('intersectValue'), width = 12
        )
      )
    )
  )# sidebarLayout
))



# ------ SERVER ------
server <- shinyServer(function(input, output) {

  data <- reactive({
    req(input$variant)
    l <- list()
    for (s in input$variant){
      l[s] = list(c(variants[[s]][["JS"]]))
    }
    as.list(l)
  })
  
  plt <- function(x){
    # --- plot the Venn diagram with the source code from package ggVennDiagram ---
    # Github: https://github.com/gaospecial/ggVennDiagram
    # Reference: https://cran.rstudio.com/web/packages/ggVennDiagram/vignettes/fully-customed.html
    venn <- Venn(x)
    venn_data <- process_data(venn)
    items <- venn_region(venn_data) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, collapse = ', '), width = 40)) %>%
      sf::st_as_sf()
    label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
    ggplot(items) +
      geom_sf(aes_string(fill="count")) +
      geom_sf_text(aes_string(label = "name"),
                   data = venn_data@setLabel,
                   inherit.aes = F) +
      geom_text(aes_string(label = "count", text = "text"),
                x = label_coord[,1],
                y = label_coord[,2],
                show.legend = FALSE) +
      theme_void() +
      scale_fill_distiller(palette = "Pastel1") + 
      scale_x_continuous(expand = expansion(mult = .2)) +
      theme(text=element_text(size=17,family="Helvetica"))
  }# ggplot object
  
  table <- function(x){
    venn <- Venn(x)
    venn_data <- process_data(venn)
    dt <- data.frame(venn_data@region, stringsAsFactors = FALSE)[c('name', 'count', 'item')]
    colnames(dt) <- c('Variant sets','Intersecting sgRNA counts','Intersecting sgRNA')
    format <- function(y){ 
      paste0(strsplit(toString(y), split = '\\.\\.')[[1]], collapse = ', ')
    }
    format2 <- function(z){ 
      paste0(strsplit(toString(unlist(z)), split = '\\,\\ ')[[1]], collapse = ', ')
    }
    dt['Variant sets'] <- apply(dt['Variant sets'], 1, format)
    dt['Intersecting sgRNA'] <- apply(dt['Intersecting sgRNA'], 1, format2)
    return(dt)
  }
  
  # --- output variables ---
  output$venn <- renderPlotly({
    # This is to show the intersect values when viewing with browser
    ax <- list(showline = FALSE)
    plotly::ggplotly(plt(data()), tooltip = c("text"), height = 800, width = 1000) %>%
      plotly::layout(xaxis = ax, yaxis = ax) # plotly object
  })
  
  output$intersectValue <- DT::renderDataTable(
    table(data()), rownames = FALSE, 
                    options = list(scrollX = TRUE,  
                                   order = list(list(1, 'desc')), 
                                   autoWidth = TRUE,
                                   columnDefs = list(list(width = '150px', targets = 0),
                                                list(width = '120px', targets = 1),
                                                list(className = 'dt-center', 
                                                     targets = 1))
                                   )
  )

  output$download1 <- downloadHandler(
    filename = function() {
      paste0('Venn_diagram_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, plot = plt(data()), device = "pdf", 
             dpi = 400, width = 17.4, height = 11.5, units = "cm")
    }
  )
  
  output$download2 <- downloadHandler(
    filename = function() {
      paste0('sgRNA_intersection_list','.csv', sep='')
    },
    content = function(file) {
      write.csv(table(data()), file)
    }
  )
  
})


#options(shiny.host = '0.0.0.0')
#options(shiny.port = 8888)
shinyApp(ui = ui, server = server)


