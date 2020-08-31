library(scmamp)
library(ggplot2)
library(scales)
library(ggrepel)
library(PMCMR)

DSAA_loadRDS_allQuantifiers <- function(){
  lf <- list.files("./result_analysis/shiny/shinyVariables/all_quantifiers/")
  all <- NULL
  nn <- NULL
  for(i in lf){
    x <- readRDS(paste0("./result_analysis/shiny/shinyVariables/all_quantifiers/",i))
    all <- c(all, list(x))
    nn <- c(nn, strsplit(strsplit(i,"_")[[1]][2],".rds")[[1]])
  }
  names(all) <- nn
  saveRDS(all, "./result_analysis/shiny/shinyVariables/all_quantifiers.rds")
  return(1)
  
}

DSAA_laod_rds <- function(){
  all <- readRDS("./result_analysis/shiny/shinyVariables/all_quantifiers.rds")
  quantifiers <- names(all)
  re <- NULL
  for(i in 1:length(quantifiers)){
    aux <- all[[i]]
    aux$Test_Size <- as.numeric(as.vector(aux$Test_Size))
    aux$R_1 <- as.numeric(as.vector(aux$R_1))
    re <- rbind(re, aux)
  }
  return(re)
}

DSAA_load_results <- function(){
  re <- DSAA_laod_rds()
  nn <- names(re)
  nn[6] <- "batchSize"
  nn[7] <- "counter"
  nn[5] <- "MAE"
  names(re) <- nn
  re$MAE <- as.numeric(as.vector(re$MAE))
  re$counter <- as.factor(re$counter)
  return(re)
}


DSAA_join_rds_results_dts_in_file <- function(counter, path_d = "results"){
  dts <- list.files("./datasets/")
  resultsTable <- NULL
  for(i in 1:length(dts)){
    qfile   <- paste0("./result_analysis/results/",path_d,"/QE_results_",counter,"_",dts[i],".rds")
    ifelse(file.exists(qfile), x <- readRDS(qfile), print(paste("NOT FOUND -- ", qfile,i,sep="----")))
    x$Qnt <- counter
    resultsTable <- rbind(resultsTable, x)
  }
  saveRDS(resultsTable, paste0("./result_analysis/shiny/shinyVariables/all_quantifiers/results_", counter, ".rds"))
}

# Run this function for including a new algorithm after running our experimental design
DSAA_include_new_algorithm <- function(new_alg_name){
  DSAA_join_rds_results_dts_in_file(new_alg_name, path_d = paste0("results_", new_alg_name))
  DSAA_loadRDS_allQuantifiers()
  re <- DSAA_load_results()
  saveRDS(re, "./result_analysis/shiny/shinyVariables/boxplot.rds")
}

# Plot figure 4
DSAA_comparisons_mixture_models <- function(){
  
  x <- readRDS("./result_analysis/shiny/shinyVariables/boxPlot.rds")
  x <- x[which(x$counter%in%c("SMM", "HDy","SORD","DyS-TS")),]
  aux <- as.data.frame(aggregate(as.numeric(x$MAE), by=list(x$counter, x$dataset), FUN=mean))
  names(aux) <- c("counter", "dataset", "MAE")
  qnt <- unique(aux$dataset)
  re <- NULL
  for(i in qnt)
    re <- rbind(re, aux[aux$dataset==i,"MAE"])
  
  re <- round(re,3)
  re <- as.data.frame(re)
  names(re) <- as.vector(aux[aux$dataset==i,"counter"])
  
  par(mar=c(0,0,0,0))
  plotCD(1-re, cex = 2.2)
  rownames(re) <- qnt
  print(nrow(re))
  re <- re[c(4,3,5,6,8,7,9,10,11,12,13,1,22,14,15,16,17,18,19,20,21,2,23,25,24),]
  return(re)
}

# Plot figure 5
DSAA_time_cost <- function(){
  vf <- list.files("./performance/results/")
  re <- NULL
  for(i in vf){
    print(paste0("./performance/results/", i))
   aux <- readRDS(paste0("./performance/results/", i))
   #aux <- aux[-which(aux$counter%in%as.vector(exc$counter)), ]
   re <- rbind(re, aux)
  }
  
  re$time <- as.numeric(as.vector(re$time))
  re$counter <- as.factor(re$counter)
  re$size <- as.numeric(as.vector(re$size))
  re$time <- re$time/1000000000
  reAll <- re
  
  # Uncomment next line for including small test set sizes 
  #re <- re[-which(re$size < 10000),]
  
  re$time <- 1/re$time
 re <- re[re$counter%in%c("DyS-TS", "SORD", "SMM", "HDy"),]
  vcol <- c("dodgerblue2",
            "indianred",
             "#A58AFF",
            "lightsteelblue",
            "plum3",
            "pink3",
            "lightsalmon2",
            "navajowhite2",
            "lightyellow4",
            "olivedrab3",
            "green3","#597DBE", "#76BF72", "#AF85BE", "#C76E6E", "#dfdbe0", "#9a989b", "#a0e09d", "#b4dded", "#078fc4", "#abdef2", "#f1b5aa"
  )
  
  
  re <- as.data.frame(re)
  p2 <- ggplot(re, aes(x=size, y=time, linetype=counter)) + 
    geom_line(aes(color=counter),size=4) +
    scale_y_continuous(trans='log2',  limits = c(0.001, 1500^2), breaks=c(0.001,1,1000, 1000^2), label = c("1m","1","1k", "1M") )+
    
    scale_linetype_manual(values=c("solid", "dotted", "dashed", "twodash"))+
    xlim(10000,100000)+
    scale_color_manual(values=vcol)+
  
    xlab("Test set size") +
    labs(title="",y="Samples per \n second (log-scale)")+
    theme_bw()+
    theme(
      #axis.text.x = element_text(angle = 0, hjust = 0),
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      #panel.grid.minor = element_blank()
      ,panel.border = element_blank()
      ,text = element_text(size=60),
      strip.background=element_rect(fill="black"),
      #legend.position = c(0.8, 0.6),
      legend.position="top",
      legend.background = element_blank(),
      legend.title = element_blank(),
      #axis.line.y = element_blank(),
      axis.line = element_line(size = 1),
      axis.text.y = element_text(size = 35),
      axis.text.x = element_text(size = 35),
      #axis.ticks.y = element_blank(),
      axis.line.x = element_line(size=1),
      legend.key.width = unit(4.5,"cm"),
      axis.ticks = element_line(size = 1),
      #,plot.margin = margin(2,2,2,2, "cm")
      #legend.text = element_text(face = "italic")
      legend.text = element_text(size = 35)
    )
  
  print(p2)
  return(reAll)
}

# Plot Figure 6
DSAA_boxPlot_rank <- function(){
  
  data <- readRDS("./result_analysis/shiny/shinyVariables/boxPlot.rds")
  data$MAE <- as.numeric(as.vector(data$MAE))
  data <- data[data$batchSize%in%c(seq(10,100,10), seq(200,500,100)),]
  df_mean <- as.data.frame(aggregate(as.numeric(data$MAE), by=list(data$counter, data$dataset), FUN=mean))
  
  names(df_mean) <- c("Counter", "dataset", "MAE")
  
  dts <- unique(df_mean$dataset)
  df_mean$MAE <- round(df_mean$MAE,3)
  df_rank <- NULL
  for(i in dts){
    x <- df_mean[df_mean$dataset==i,]
    aux <- as.data.frame(cbind(as.vector(x[,1]), rank(as.numeric(as.vector(x[,3])))))
    aux$dataset <- i
    df_rank <- rbind(df_rank, aux)
  }
  
  df_rank <- as.data.frame(df_rank)
  names(df_rank) <- c("Counter", "Rank", "Dataset")
  
  df_rank$Rank <- as.numeric(as.vector(df_rank$Rank))
  print(length(unique(df_rank$Dataset)))
  theme_set(theme_bw())
  v1 <- ggplot(df_rank, aes(reorder(Counter, Rank), Rank, fill=Counter)) + 
    geom_boxplot(varwidth=F, alpha=0.2, lwd=1.2)+
    labs(title="",x="Quantification algorithms")+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=35),
      plot.background = element_blank()
      #,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank()
      ,text = element_text(size=50),
      strip.background=element_rect(fill="black"),
      #legend.position = c(.9, 0.9),
      legend.position="none",
      #legend.position="top",
      legend.background = element_blank()
      ,legend.title = element_blank()
      ,legend.key.width = unit(1,"cm"),
      
      #legend.background = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(size=1)
      #legend.key.width = unit(3.5,"cm"),
      #axis.ticks = element_line(size = 1)
      #,plot.margin = margin(2,2,2,2, "cm")
      #legend.text = element_text(face = "italic")
    )
  
  print(v1)
  
  return(df_rank)
  
}


# plot Figure 1
DSAA_intrXrank <- function(){
  re <- DSAA_boxPlot_rank()
  rt <- DSAA_time_cost()
  
  rt$time <- rt[,3]/rt[,2]
  rt <- rt[,-3]
  aux1 <- aggregate(re$Rank, by=list(re$Counter), FUN=mean)
  aux2 <- aggregate(rt$time, by=list(rt$counter), FUN=mean)
  
  names(aux1) <- c("Counter", "Rank")
  names(aux2) <- c("Counter", "Speed")
  
  #rt <- rt[-which(rt$size%in%c(100,1000)),]
  
  if(sum(aux1$Counter=="QT")){
    pos <- which(aux1$Counter=="QT")
    aux2 <- as.matrix(aux2)
    aux2 <- rbind(aux2[1:(pos-1),], c("QT", 0), aux2[pos:nrow(aux2),])
  }
  
  
  aux2 <- as.data.frame(aux2)
  aux2$Speed[aux2$Counter=="QT"] <- aux2$Speed[aux2$Counter=="CC"]
  
  df <- cbind(aux1,aux2[,-1])
  names(df) <- c("Counter", "Rank", "Speed")
  
  df$Speed <- as.numeric(as.vector(df$Speed))
  
  tmp <- with(df,data.frame(x=c(1, 1, 6), y=c(30000, 1100^3, 1100^3)))
  
  theme_set(theme_bw())
  v1 <- ggplot(df) + 
    geom_point(aes(x=Rank, y = Speed, colour=Counter), size=3)+
    
    scale_y_continuous(trans='log2',  limits = c(500, 1500^3), breaks=c(1000, 1000^2, 1000^3), label = c("1k", "1M", "1G") )+
    
    geom_text_repel(aes(Rank, Speed, label = Counter, colour=Counter), size = 13)+
    geom_polygon(data=tmp, aes(x=x, y=y),fill="#d8161688")+
    labs(title="",y="Instances per \n second (log-scale)")+
    
    scale_x_continuous( limits=c(1, 14), breaks=seq(1,14,2)) +
    
    # quadrant labels
    annotate("text", x = 1.5, y = 26.5, alpha = 0.35, label = "Fast and", size=8) +
    annotate("text", x = 1.5, y = 10.1, alpha = 0.35, label = "accurate", size=8) +
    
    annotate("text", x = 1.5, y = 0.0005, alpha = 0.35, label = "Slow and", size=8) +
    annotate("text", x = 1.5, y = 0.00018, alpha = 0.35, label = "accurate", size=8) +
    
    annotate("text", x = 13.5, y = 0.0005, alpha = 0.35, label = "Slow and", size=8) +
    annotate("text", x = 13.5, y = 0.00018, alpha = 0.35, label = "inaccurate", size=8) +
    
    annotate("text", x = 13.5, y = 26.5, alpha = 0.35, label = "Fast and", size=8) +
    annotate("text", x = 13.5, y = 10.1, alpha = 0.35, label = "inaccurate", size=8) +
    
    theme(
      #axis.text.x = element_text(angle = 0, size = 25),
      #axis.text.y = element_text(angle = 90, size=25),
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank()
      ,text = element_text(size=50),
      strip.background=element_rect(fill="black"),
      #legend.position = c(.9, 0.9),
      legend.position="none",
      #legend.position="top",
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.box = "horizontal",
      #,legend.key.width = unit(1,"cm"),
      
      #legend.background = element_blank(),
      #axis.text.y = element_blank(),
      #axis.ticks.y = element_blank(),
      #axis.line.y = element_blank(),
      axis.line = element_line(size = 1),
      
      #axis.line.x = element_blank(),
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank()
      
      #legend.key.width = unit(3.5,"cm"),
      axis.ticks = element_line(size = 1)
      #,plot.margin = margin(2,2,2,2, "cm")
      #legend.text = element_text(face = "italic")
    )
  print(v1)
  
  #ggsave(file="runtimeXrank.svg", plot=v1, width=10, height=13)
  
  print(min(df$Speed))
  
  return(df)
}

# plot Figure 7
DSAA_heapmap <- function(df_heat){
  
  if(is.null(df_heat)){
    vcol <- c("dodgerblue2",
              "indianred",
              "#A58AFF",
              "lightsteelblue",
              "plum3",
              "pink3",
              "lightsalmon2",
              "navajowhite2",
              "lightyellow4",
              "olivedrab3",
              "green3",
               "#b4dded", 
              "#078fc4",  
              "darkgoldenrod3",
              "#abdef2",
              "#76BF72",
              "#597DBE",
              "#AF85BE",
              "#C76E6E",
              "#dfdbe0", 
              "#9a989b", 
              "#a0e09d", 
              "#f1b5aa"
    )
    
    x <- readRDS("./result_analysis/shiny/shinyVariables/boxPlot.rds")
    re<- x[x$batchSize%in%c(seq(10,100,10),seq(200,500,100)),]
    df_heat <- aggregate(as.numeric(re$MAE), by=list(re$counter, re$batchSize), FUN=mean)
    df_heat <- as.data.frame(df_heat)
    names(df_heat) <- c("counter", "size", "error")
    
    df_heat$size <- as.factor(df_heat$size)
    siz <- unique(df_heat$size)
    re <- NULL
    co <- c("DyS-TS", "SORD", "SMM", "MAX", "EMQ", "ACC", "X", "MS2", "MS", "PACC", "CC", "PCC","HDy", "T50", "QT")
    for(i in siz){
      aux <- df_heat[df_heat$size==i,]
      for(j in co) re <- rbind(re, aux[aux$counter==j,])
    }
    df_heat <- re
    df_heat$counter <- factor(df_heat$counter, levels=co)
  }
  
  co <- c("DyS-TS", "SORD", "SMM", "MAX", "EMQ", "ACC", "X", "MS2", "MS", "PACC", "CC", "PCC","HDy", "T50", "QT")
  colfunc<-colorRampPalette(c("#b50600", "#f0ec07","lightsalmon2","green3","#0015cf"))
  
  df_heat$error <- as.numeric(as.vector(df_heat$error))
  
  theme_set(theme_bw())
  
  jco <-  c("default", "pal2", "pal3", "pal4", "pal5",
            "pal6", "pal7", "pal8", "pal9", "pal10", "pal11", "pal12", "rainbow")
  
  p<- ggplot(df_heat, aes(y=error,x=counter)) +
    scale_size_continuous(range = c(4, 18), guide = F)+
    geom_point(aes(fill=size), pch=21, stroke=1, alpha=1, size=10)+
    scale_fill_manual(values =colfunc(15))+
    scale_y_continuous(trans='log2')+
    labs(fill = "")+
    ylab("MAE \n (log-scale)") +
    labs(title="",x="Quantification algorithms")+
    theme(
      plot.background = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank()
      ,text = element_text(size=50)
      ,legend.position = "top"
      ,axis.ticks.y = element_blank()
      ,axis.line.x = element_line(size=1.5)
      ,axis.ticks = element_line(size = 2)
      ,axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      
    )+guides(fill = guide_legend(nrow = 2), size=F) #+guides(fill = guide_legend(title.position = "top"))
    print(p)
    return(df_heat) 
}

# Table IV
DSAA_stats_comparisons <- function(){
  
  x <- readRDS("./result_analysis/shiny/shinyVariables/boxPlot.rds")
  x <- x[x$batchSize%in%c(seq(10,100,10),seq(200,500,100)),]
  
  aux <- as.data.frame(aggregate(as.numeric(x$MAE), by=list(x$counter, x$dataset, x$batchSize), FUN=mean))
  names(aux) <- c("counter", "dataset", "size","MAE")
  
  result <- NULL
  
  dts <- unique(aux$dataset)[c(4,3,5,6,8,7,9,10,11,12,13,1,22,14,15,16,17,18,19,20,21,2,23,25,24)]
  qnt <- c("DyS-TS", "SORD", "SMM", "MAX", "EMQ", "ACC", "X", "MS2", "MS", "PACC", "CC", "PCC","HDy", "T50", "QT")
  for(i in dts){
    re <- NULL
    print(i)
    for(j in qnt){
      re <- cbind(re, aux[aux$dataset==i & aux$counter==j,"MAE"])
    }
    
    re <- as.data.frame(re)
    names(re) <- qnt
    print(posthoc.friedman.nemenyi.test(apply(re, 2,as.numeric)))
    print("####################################################################")
    result <- rbind(result, apply(re, 2, mean))
  }
  return(result)
  
}

#-------------------------------------------------------------------------------------
# Figure 1
DSAA_intrXrank()

# Figure 4
DSAA_comparisons_mixture_models()

# Figure 5
DSAA_time_cost()

# Figure 6
DSAA_boxPlot_rank()

# Figure 7
DSAA_heapmap()

#Table IV and statistical comparisons are printed 
tableIV <- round(DSAA_stats_comparisons(),3)


