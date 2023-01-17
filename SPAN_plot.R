library(optparse)
## Check arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="file name", metavar="character"),
  make_option(c("-f", "--main"), type="character", default=NULL, 
              help="main file", metavar="character"),
  make_option(c("-n", "--number_of_samples"), type="numeric", default=4, 
              help="number_of_samples filter", metavar="numeric"),
  make_option(c("-r", "--rpm"), type="numeric", default=5000, 
              help="reads per milliom filter", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied", call.=FALSE)
}
## Import library and working path
library(tidyverse)
library(ggplot2)
library(ftplottools)
library(ggpubr)
library(cowplot)
path <- getwd()
file_path <- paste(path,"/raw_tsv/",sep= "")
## The ncsgRNAs distribution plot
## ggplot genome
rect_df <- data.frame(xmin = c(0,266,13483,21562,25392,26244,26522,27201,27393,27755,27893,28273,29558),
                      xmax = c(265,13483,21555,25384,26220,26472,27191,27387,27759,27887,28259,29533,
                               29674),
                      names =factor(
                        c('5-UTR','ORF1a','ORF1b','S','3a','E','M','6','7a','7b','8','N','3-UTR')))
Genome_annotation = ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.02, ymax = 0.02, fill = names), 
            alpha = 1,data = rect_df,inherit.aes = FALSE) +
  scale_fill_manual(name='ORF',
                    values=c('5-TUR'='#E3E3E3',
                             'ORF1a'='#EADDBD', 
                             'ORF1b'='#E8E1CF', 
                             'S'='#FCF8C1',
                             '3a'='#DEFCCE',
                             'E'='#E9DBF9',
                             'M'='#CEEEF9',
                             '6'='#FCE2C1',
                             '7a'='#FCCEEB',
                             '7b'='#DDDBF9',
                             '8'='#BCD8E1',
                             'N'='#E1BCBC',
                             '3-UTR'='#E3E3E3'),
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
## Import Data and Summarize
Main = read.delim(file.path(paste(file_path,opt$main,"_sgRNA_nc_junction_summary.tsv",sep ="")),
           sep="\t",
           col.names =c("j5","j3","count","startpos",
                        "productsize","JS","deletion size",
                        "sample","sgRPM","chrom","type","name",
                        "start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'))

Main_junction <-Main %>%
  group_by(JS,j5,j3)%>%
  summarise(number_of_reads = sum(count),RPM=sum(sgRPM),number_of_samples=n(),
            variants=opt$main,y=0.02)
sample_list <- read.delim(file.path(paste(path,opt$input,sep ="/")), sep="\n",header = TRUE)
for (i in sample_list$Lineages){
  Lineage = read.delim(file.path(paste(file_path,i,"_sgRNA_nc_junction_summary.tsv",sep ="")),
                         sep="\t",
                         col.names =c("j5","j3","count","startpos",
                                      "productsize","JS","deletion size",
                                      "sample","sgRPM","chrom","type","name",
                                      "start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'))
  Lineage_junction <- Lineage %>%
    group_by(JS,j5,j3)%>%
    summarise(number_of_reads = sum(count),RPM=sum(sgRPM),number_of_samples=n(),
              variants=i,y=-0.02) 
  Variant = rbind(Main_junction,Lineage_junction)
  Fig=Genome_annotation+
    geom_curve(data=Variant%>%filter(variants==i)%>%filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm),lineend = "round",curvature=0.2,size=0.1,color="#464646")+
    geom_curve(data=Variant%>%filter(variants==opt$main)%>%filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm),lineend = "round",curvature=-0.2,size=0.1,color="#464646")+
    geom_point(alpha=0.8,data=Variant%>%filter(variants==i)%>%filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%arrange(number_of_samples),aes(x=j5,y=-0.02,size=RPM,color=number_of_samples))+
    geom_point(alpha=0.8,data=Variant%>%filter(variants==opt$main)%>%filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%arrange(number_of_samples),aes(x=j5,y=0.02,size=RPM,color=number_of_samples))+
    geom_point(alpha=0.8,data=Variant%>%filter(variants==i)%>%filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%arrange(number_of_samples),aes(x=j3,y=-0.02,size=RPM,color=number_of_samples))+
    geom_point(alpha=0.8,data=Variant%>%filter(variants==opt$main)%>%filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%arrange(number_of_samples),aes(x=j3,y=0.02,size=RPM,color=number_of_samples))+
    scale_color_gradient("Number of\nSamples",low="#E3D7D5",high="#EA5F51")+
    annotate("text",label=i,fontface='bold',x=1000,y=-0.05,color="#000000", hjust = 0,family = "Helvetica")+
    annotate("text",label=opt$main,fontface='bold',x=1000,y=0.05,color="#000000", hjust = 0,family = "Helvetica")+
    scale_size_continuous("RPM",limits= c(0,5000))+
    theme(legend.position = 'none',
          legend.title = element_blank())
  rect_plot <- ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y))+
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.02, ymax = 0.02, fill = names), 
              alpha = 1,data = rect_df,inherit.aes = FALSE) +
    scale_fill_manual(name='ORF',
                      values=c('5-TUR'='#E3E3E3',
                               'ORF1a'='#EADDBD', 
                               'ORF1b'='#E8E1CF', 
                               'S'='#FCF8C1',
                               '3a'='#DEFCCE',
                               'E'='#E9DBF9',
                               'M'='#CEEEF9',
                               '6'='#FCE2C1',
                               '7a'='#FCCEEB',
                               '7b'='#DDDBF9',
                               '8'='#BCD8E1',
                               'N'='#E1BCBC',
                               '3-UTR'='#E3E3E3'),
                      breaks = c('ORF1a','ORF1b','S','3a','E','M','6','7a','7b','8','N'))+
    theme(legend.position = 'top',
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.height = unit(0.1, "cm"),
          legend.key.width = unit(0.5, "cm"))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  rect_legend <- ggpubr::get_legend(rect_plot)
  
  color_plot <-ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y)) + 
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==i)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j5,y=-0.02,color=number_of_samples))+
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==opt$main)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j5,y=0.02,color=number_of_samples))+
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==i)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j3,y=-0.02,color=number_of_samples))+
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==opt$main)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j3,y=0.02,color=number_of_samples))+
    scale_color_gradient("Number of\nSamples",low="#E3D7D5",high="#EA5F51",limits= c(0,20))+
    theme(legend.position = "right",
          legend.direction = "vertical",
          legend.key.height = unit(0.4, "cm"))
  color_legend <- ggpubr::get_legend(color_plot)
  
  size_plot <-ggplot(mapping = aes(x=j5,xend=j3,y=y,yend=y)) + 
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==i)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j5,y=-0.02,size=RPM))+
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==opt$main)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j5,y=0.02,size=RPM))+
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==i)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j3,y=-0.02,size=RPM))+
    geom_point(alpha=0.8,data=Variant%>%
                 filter(variants==opt$main)%>%
                 filter(number_of_samples>opt$number_of_samples,RPM>opt$rpm)%>%
                 arrange(number_of_samples),aes(x=j3,y=0.02,size=RPM))+
    scale_size_continuous("RPM",limits= c(0,5000))+
    theme(legend.position = "right",
          legend.direction = "vertical",
          legend.key.height = unit(0.4, "cm"))
  size_legend <- ggpubr::get_legend(size_plot)
  # combine all these elements
  Fig<-cowplot::plot_grid(plotlist = list(rect_legend,NULL, NULL,Fig, color_legend,size_legend),
                            rel_heights = c(1, 5),
                            rel_widths =  c(4, 1, 1)) 
  ggsave(file.path(
    paste(i,opt$main,"ncsgRNAs distribution plot among different Omicron lineages.pdf",sep ="_")), 
                   Fig, width = 17.4, height = 5.75, units = "cm")
}

## For loop automatically generate all intersecting ncsgRNAs plots (separately)
## Import Data
All_Omicron_intersect = read.delim(file.path(paste(path,"/All_intersect_Omicron_ncsg.tsv",sep ="")),
                  sep="\t",
                  col.names =c("lineage","j5","j3","count","startpos",
                "productsize","JS","deletion size",
                "sample","sgRPM","chrom","type","name",
                "start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'))
## Main function
All_Omicron_intersect <- All_Omicron_intersect %>%group_by(JS)
All_key <- All_Omicron_intersect %>% group_keys(JS)
for (i in All_key$JS){
  plot_sg<-All_Omicron_intersect %>% filter(JS == i)
  plot<-ggplot(plot_sg, aes(x = as.factor(lineage), y = log10(sgRPM),color=as.factor(lineage))) +
  geom_boxplot() + 
  geom_jitter(shape = 16,                        
              position = position_jitter(0.2)) + 
  xlab("") +   
  ylab(paste(i," RPM(log10)")) +   
  labs(color = "lineage") +
  coord_cartesian(ylim = c(-1, 3))+
  scale_fill_distiller(palette = "Pastel1")+
  #scale_color_manual(values = c("B.1.1.529"="#ffe6e6", 
  #                              "BA.1"="#fcf8c1",
  #                              "BA.2"="#defcce",
  #                              "BA.4"="#e9dbf9",
  #                              "BA.5"="#ceeef9"))+
  ft_theme()+
  theme(legend.position="none")+
  stat_compare_means(label.y = 3)
  #stat_compare_means(label.y = 3,label.x = "BA.1") 
  ggsave(file.path(paste(i,"ncsgRNAs expression among different Omicron lineages.pdf")), plot, width = 17.4, height = 23, units = "cm",dpi=700)
}