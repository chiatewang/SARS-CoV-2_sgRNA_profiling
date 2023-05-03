library(shiny)
library(rjson)
library(ggVennDiagram)
library(ggplot2)
library(plotly)
library(tidyverse)
library(ftplottools)
library(ggpubr)
library(cowplot)
library(DT)
path <- getwd()
# ------ Read in the JSON file to get Junction Site ------
variants <- fromJSON(file = file.path(paste(path,"/sgRNA_list.json",sep = "")))
# Ref. https://zhuanlan.zhihu.com/p/148779322
# ------ UI ------
ui <- shinyUI(fluidPage(
  titlePanel("Profiling noncanonical subgenomic RNAs of SARS-CoV-2"),
  tabsetPanel(
    tabPanel("Venn Diagram",
             sidebarPanel(width = 3,
                          checkboxGroupInput("variant", label = "Variant:",
                                             choiceNames = names(variants),
                                             choiceValues = names(variants)),
                          downloadButton("download1", "Download Venn Diagram Graph"),
                          downloadButton("download2", "Download Intersecting Table")),
             mainPanel(width = 9,
                       fluidRow(plotlyOutput("venn", width = "100%"), style = "height:800px"),
                       fluidRow(column(dataTableOutput('intersectValue'), width = 12))
             )
    ),
    tabPanel("Junction Site Sashimi Plot",
             sidebarPanel(width = 3,
                          selectInput("variant1", label = "Variant1:",
                                      choice = names(variants)),
                          selectInput("variant2", label = "Variant2:",
                                      choice = names(variants)),
                          sliderInput(inputId = "RPM",
                                      label = "RPM >=",
                                      min = 0,
                                      max = 50,
                                      value = 1),
                          sliderInput(inputId = "ns",
                                      label = "Number of samples >=",
                                      min = 0,
                                      max = 10,
                                      value = 0),
                          downloadButton("download3", "Download Junction Sashimi Plot")),
             mainPanel(width = 9,
                       plotOutput("junction"),
             )
    ),
    tabPanel("Boxplot",
             sidebarPanel(width = 3,
                          fileInput("file", 
                                    label = "Please upload intersecting CSV: ", 
                                    accept = ".csv"),
                          #checkboxInput("header", "Header", TRUE),
                          #downloadButton("download4", "Download Boxplot")
                          ),
             mainPanel(width = 9,
                       uiOutput("boxplot")
             )
    ))
  )
)

# ------ Server ------
server <- shinyServer(function(input, output){
  # ------ Function 1 Venn Diagram ------
  data_JS <- reactive({
    req(input$variant)
    l_JS <- list()
    for (js in input$variant){
      l_JS[js] = list(c(variants[[js]][["JS"]]))
    }
    as.list(l_JS)
  })
  venndiagram <- function(x){
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
  }
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
  # ------ Function 2 Junction Sashimi Plot ------
  data_JS1 <- reactive({
    req(input$variant1)
    l_JS1 <- list()
    for (js in input$variant1){
      l_JS1[js] = list(c(variants[[js]][["JS"]]))
    }
    l_JS1 = unlist(l_JS1)
  })
  data_nr1 <- reactive({
    req(input$variant1)
    l_nr1 <- list()
    for (nr in input$variant1){
      l_nr1[nr] = list(c(variants[[nr]][["number_of_reads"]]))
    }
    l_nr1 = unlist(l_nr1)
  })
  data_RPM1 <- reactive({
    req(input$variant1)
    l_RPM1 <- list()
    for (RPM in input$variant1){
      l_RPM1[RPM] = list(c(variants[[RPM]][["RPM"]]))
    }
    l_RPM1 = unlist(l_RPM1)
  })
  data_ns1 <- reactive({
    req(input$variant1)
    l_ns1 <- list()
    for (ns in input$variant1){
      l_ns1[ns] = list(c(variants[[ns]][["number_of_samples"]]))
    }
    l_ns1 = unlist(l_ns1)
  })
  data_j51 <- reactive({
    req(input$variant1)
    l_j51 <- list()
    for (j5 in input$variant1){
      l_j51[j5] = list(c(variants[[j5]][["j5"]]))
    }
    l_j51 = unlist(l_j51)
  })
  data_j31 <- reactive({
    req(input$variant1)
    l_j31 <- list()
    for (j3 in input$variant1){
      l_j31[j3] = list(c(variants[[j3]][["j3"]]))
    }
    l_j31 = unlist(l_j31)
  })
  
  data_JS2 <- reactive({
    req(input$variant2)
    l_JS2 <- list()
    for (js in input$variant2){
      l_JS2[js] = list(c(variants[[js]][["JS"]]))
    }
    l_JS2 = unlist(l_JS2)
  })
  data_nr2 <- reactive({
    req(input$variant2)
    l_nr2 <- list()
    for (nr in input$variant2){
      l_nr2[nr] = list(c(variants[[nr]][["number_of_reads"]]))
    }
    l_nr2 = unlist(l_nr2)
  })
  data_RPM2 <- reactive({
    req(input$variant2)
    l_RPM2 <- list()
    for (RPM in input$variant2){
      l_RPM2[RPM] = list(c(variants[[RPM]][["RPM"]]))
    }
    l_RPM2 = unlist(l_RPM2)
  })
  data_ns2 <- reactive({
    req(input$variant2)
    l_ns2 <- list()
    for (ns in input$variant2){
      l_ns2[ns] = list(c(variants[[ns]][["number_of_samples"]]))
    }
    l_ns2 = unlist(l_ns2)
  })
  data_j52 <- reactive({
    req(input$variant2)
    l_j52 <- list()
    for (j5 in input$variant2){
      l_j52[j5] = list(c(variants[[j5]][["j5"]]))
    }
    l_j52 = unlist(l_j52)
  })
  data_j32 <- reactive({
    req(input$variant2)
    l_j32 <- list()
    for (j3 in input$variant2){
      l_j32[j3] = list(c(variants[[j3]][["j3"]]))
    }
    l_j32 = unlist(l_j32)
  })
  Variant <- reactive({
    Variant1<-do.call(data.frame,list('variants' = input$variant1,
                                      'JS' = data_JS1(),
                                      'number_of_reads' = data_nr1(),
                                      'RPM' = data_RPM1(),
                                      'number_of_samples' = data_ns1(),
                                      'j5' = data_j51(),
                                      'j3'= data_j31(),
                                      'y'=0.02))
    Variant2<-do.call(data.frame,list('variants' = input$variant2,
                                      'JS' = data_JS2(),
                                      'number_of_reads' = data_nr2(),
                                      'RPM' = data_RPM2(),
                                      'number_of_samples' = data_ns2(),
                                      'j5' = data_j52(),
                                      'j3' = data_j32(),
                                      'y'=-0.02))
    Variant<- rbind(Variant1,Variant2)
  })
  junction <- function(){
    rect_df <- data.frame(xmin = c(0,266,13483,21562,25392,26244,26522,
                                   27201,27393,27755,27893,28273,29558),
                          xmax = c(265,13483,21555,25384,26220,26472,27191,
                                   27387,27759,27887,28259,29533,29674),
                          names =factor(c('5-UTR','ORF1a','ORF1b','S','3a','E',
                                          'M','6','7a','7b','8','N','3-UTR')))
    Genome_annotation = ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y))+
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.02, ymax = 0.02, fill = names), 
                alpha = 1,data = rect_df,inherit.aes = FALSE) +
      scale_fill_manual(name='ORF',
                        values=c('5-TUR'='#E3E3E3','ORF1a'='#EADDBD', 'ORF1b'='#E8E1CF', 
                                 'S'='#FCF8C1','3a'='#DEFCCE','E'='#E9DBF9','M'='#CEEEF9',
                                 '6'='#FCE2C1','7a'='#FCCEEB','7b'='#DDDBF9','8'='#BCD8E1',
                                 'N'='#E1BCBC','3-UTR'='#E3E3E3'),
                        breaks = c('ORF1a','ORF1b','S','3a','E','M','6','7a','7b','8','N'))+
      scale_x_continuous("Position")+
      scale_alpha_continuous("Number of Samples")+
      ft_theme()+
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.line.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.direction = "vertical", 
            legend.box = "vertical",
            legend.position='right',
            legend.justification = "center",
            legend.key.height = unit(1, "cm"),
            legend.key.width = unit(0.2, "cm"))
    Fig<-Genome_annotation+
      geom_curve(data=Variant()%>%filter(variants==input$variant1)
                 %>%filter(number_of_samples>=input$ns,RPM>input$RPM),
                 lineend = "round",curvature=-0.1,size=0.1,color="#464646")+
      geom_curve(data=Variant()%>%filter(variants==input$variant2)
                 %>%filter(number_of_samples>=input$ns,RPM>input$RPM),
                 lineend = "round",curvature=0.1,size=0.1,color="#464646")+
      geom_point(alpha=0.8,data=Variant()%>%filter(variants==input$variant1)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j5,y=0.02,size=RPM,color=number_of_samples))+
      geom_point(alpha=0.8,data=Variant()%>%filter(variants==input$variant2)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j5,y=-0.02,size=RPM,color=number_of_samples))+
      geom_point(alpha=0.8,data=Variant()%>%filter(variants==input$variant1)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j3,y=0.02,size=RPM,color=number_of_samples))+
      geom_point(alpha=0.8,data=Variant()%>%filter(variants==input$variant2)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j3,y=-0.02,size=RPM,color=number_of_samples))+
      scale_color_gradient("Number of\nSamples",low="#E3D7D5",high="#EA5F51")+
      scale_size_continuous("RPM",limits= c(0,50))+
      annotate("text",label=input$variant1,fontface='bold',x=1000,y=0.04,color="#000000")+
      annotate("text",label=input$variant2,fontface='bold',x=1000,y=-0.04,color="#000000")+
      theme(legend.position = 'none',
            legend.title = element_blank())
    rect_plot <- ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y))+
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.02, ymax = 0.02, fill = names), 
                alpha = 1,data = rect_df,inherit.aes = FALSE) +
      scale_fill_manual(name='ORF',
                        values=c('5-TUR'='#E3E3E3','ORF1a'='#EADDBD', 'ORF1b'='#E8E1CF', 
                                 'S'='#FCF8C1','3a'='#DEFCCE','E'='#E9DBF9',
                                 'M'='#CEEEF9','6'='#FCE2C1','7a'='#FCCEEB',
                                 '7b'='#DDDBF9','8'='#BCD8E1','N'='#E1BCBC','3-UTR'='#E3E3E3'),
                        breaks = c('ORF1a','ORF1b','S','3a','E','M','6','7a','7b','8','N'))+
      theme(legend.position = 'top',
            legend.direction = "horizontal",legend.box = "horizontal",
            legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.65, "cm"))+
      guides(fill=guide_legend(nrow=1,byrow=TRUE))
    rect_legend <- ggpubr::get_legend(rect_plot)
    
    color_plot <-ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y)) + 
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant1)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j5,y=0.02,color=number_of_samples))+
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant2)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j5,y=-0.02,color=number_of_samples))+
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant1)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j3,y=0.02,color=number_of_samples))+
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant2)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j3,y=-0.02,color=number_of_samples))+
      scale_color_gradient("Number of\nSamples",low="#E3D7D5",high="#EA5F51",limits= c(0,20))+
      theme(legend.position = "right",
            legend.direction = "vertical",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.15, "cm"))
    color_legend <- ggpubr::get_legend(color_plot)
    size_plot <-ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y)) + 
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant1)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j5,y=0.02,size=RPM))+
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant2)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j5,y=-0.02,size=RPM))+
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant1)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j3,y=0.02,size=RPM))+
      geom_point(alpha=0.8,data=Variant()%>%
                   filter(variants==input$variant2)%>%
                   filter(number_of_samples>=input$ns,RPM>input$RPM)%>%
                   arrange(number_of_samples),aes(x=j3,y=-0.02,size=RPM))+
      scale_size_continuous("RPM",limits= c(0,50))+
      theme(legend.position = "right",
            legend.direction = "vertical",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.15, "cm"))
    size_legend <- ggpubr::get_legend(size_plot)
    # combine all these elements
    Fig1<-cowplot::plot_grid(plotlist = list(rect_legend,NULL,NULL,Fig, size_legend,color_legend),
                             rel_heights = c(1, 10),
                             rel_widths =  c(10, 1.25, 1.25))
    plot(Fig1)
  }
  # ------ Function 3 ncsgRNA expression Boxplot ------
  summary<-reactive({
    req(input$file)
    file <- read.csv(input$file$datapath, header = TRUE)
    variant_sets = file[ nrow(file) , 2]
    intersect_filter = file[file$Variant.sets == variant_sets,4]
    intersect_filter = strsplit(intersect_filter, ", ")[[1]]
    variant_sets = strsplit(variant_sets, ", ")[[1]]
    summary<- NULL
    for (i in variant_sets){
      variant = read.table(file.path(paste(path,"/raw_tsv/",i,"_sgRNA_nc_junction_summary.tsv",sep ="")),
                           sep="\t")
      variant['lineage'] = i 
      summary=rbind(summary,variant)
    }
    colnames(summary) <- c("j5","j3","count","JS","deletion size",
                           "sample","RPM","startpos","lineage")
    return(summary)
  })
  summary2<- reactive({
    req(input$file)
    file <- read.csv(input$file$datapath, header = TRUE)
    variant_sets = file[ nrow(file) , 2]
    intersect_filter = file[file$Variant.sets == variant_sets,4]
    intersect_filter = strsplit(intersect_filter, ", ")[[1]]
    variant_sets = strsplit(variant_sets, ", ")[[1]]
    summary2<-NULL
    for (k in intersect_filter){
      summary1<-summary() %>% filter(JS == k)
      summary2<-rbind(summary2,summary1)
    }
    summary2 <- summary2 %>%group_by(JS)
    return(summary2)
  })
  All_key <- reactive({
    req(input$file)
    file <- read.csv(input$file$datapath, header = TRUE)
    variant_sets = file[ nrow(file) , 2]
    intersect_filter = file[file$Variant.sets == variant_sets,4]
    intersect_filter = strsplit(intersect_filter, ", ")[[1]]
    variant_sets = strsplit(variant_sets, ", ")[[1]]
    summary2<-NULL
    for (k in intersect_filter){
      summary1<-summary() %>% filter(JS == k)
      summary2<-rbind(summary2,summary1)
    }
    summary2 <- summary2 %>%group_by(JS)
    All_key <- summary2 %>% group_keys(JS)
    return(All_key)
  })
  Box <- function(){
    lapply(All_key()$JS, function(i) {
      output[[paste0("boxplot_", i)]] <- renderPlotly({
        plot_sg <- summary2() %>% filter(JS == i)
        plot <- ggplot(plot_sg, 
                       aes(x = as.factor(lineage), 
                           y = log10(as.numeric(RPM)),
                           color = as.factor(lineage))) +
          geom_boxplot() +
          geom_jitter(shape = 16, position = position_jitter(0.2)) +
          xlab("") +
          ylab(paste(i, " log10(RPM)")) +
          labs(color = "lineage") +
          coord_cartesian(ylim = c(-1, 3)) +
          scale_fill_distiller(palette = "Pastel1") +
          ft_theme() +
          theme(legend.position = "none") +
          stat_compare_means(label.y = 3)
        plotly::ggplotly(plot)
      })
    })
  }
  # --- output variables ---
  # ------ Output 1 Venn Diagram ------
  output$venn <- renderPlotly({
    # This is to show the intersect values when viewing with browser
    ax <- list(showline = FALSE)
    plotly::ggplotly(venndiagram(data_JS()), tooltip = c("text"), height = 800, width = 1000) %>%
      plotly::layout(xaxis = ax, yaxis = ax) 
  })
  output$intersectValue <- DT::renderDataTable(
    table(data_JS()), rownames = FALSE, 
    options = list(scrollX = TRUE,  
                   order = list(list(1, 'desc')), 
                   autoWidth = TRUE,
                   columnDefs = list(list(width = '150px', targets = 0),
                                     list(width = '120px', targets = 1),
                                     list(className = 'dt-center', 
                                          targets = 1))
    )
  )
  # ------ Output 2 Junction Sashimi Plot ------
  output$junction <- renderPlot({
    print(junction())
  })
  # ------ Output 3 ncsgRNA expression Boxplot ------
  output$boxplot <- renderUI({
    print(Box())
  })
  # ------ Download Function  ------
  output$download1 <- downloadHandler(
    filename = function() {
      paste0('Venn_diagram_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, plot = venndiagram(data_JS()), device = "pdf", 
             dpi = 400, width = 17.4, height = 5.75, units = "cm")
    }
  )
  output$download2 <- downloadHandler(
    filename = function() {
      paste0('sgRNA_intersection_list','.csv', sep='')
    },
    content = function(file) {
      write.csv(table(data_JS()), file)
    }
  )
  output$download3 <- downloadHandler(
    filename = function() {
      paste0('Junction_Sashimi_Plot_', Sys.Date(), '.pdf', sep='')
    },
    content=function(file){
      ggsave(file,plot = junction(), 
             device = "pdf",dpi = 400, width = 21.0, height = 7.425, units = "cm")
      
    }
  )

})

shinyApp(ui = ui, server = server)