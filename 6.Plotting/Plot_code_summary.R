#fig2 a
library(ggplot2)
library(ggpointdensity)
library(ggpubr)

setwd('/R_plot/fig2/a')
valid = read.csv('OCT_age_prediction_evaluation.csv') # load the data

df <- valid[, c("df.Age", "df.pre_age")]
colnames(df) <- c("CA", "BA") # choose the chronological age column and biological age column
mae <- mean(abs(df$CA-df$BA)) # calsulate the MAE
p <- ggplot(df, aes(CA,  BA)) + geom_point()+ geom_pointdensity() +
  scale_color_gradientn(colors = c('#555AA6', '#74C5A2', '#FCF8B5', '#F57948', '#AA1246'))+
  theme_classic() + geom_smooth(method = "lm", se=T,  color="black")+ 
  stat_cor(aes(label = paste(..r.label.., round(mae,2), sep = '~`,`~')), 
           method = 'pearson',  label.x.npc = 'left', label.y.npc = 'top', size = 2.7)+
  xlim(40,70)+ylim(40,70) + theme(legend.position = "none")+
  theme( axis.text.y = element_text(size=12, color="black"),  
         axis.text.x = element_text(angle = 0,hjust=1,size=12, color="black"))+ 
  labs(x="Chronological age", y="Brain_GM age", title=paste0("All participants"))

ggsave("OCT.pdf", plot = p, width = 8, height = 8) # save the plot


#fig2 b
library(ggplot2)

setwd('/R_plot/fig2/b')
data = read.csv('mae.csv')

ggplot(data, aes(x = organ, group = sex, fill = sex)) +
  geom_bar(aes(y = MAE),
           stat = "identity", 
           width = 0.5, 
           position = 'dodge' # 'dodge' means grouping, and the value 'stack' means stacking
  ) +
  scale_fill_manual(values = c('#fd3f4e',  '#0494f2')) +
  xlab("organ") +
  labs(fill = "sex")+ # color denpends on sex
  theme_minimal()


#fig2 c
library(ggplot2)
library(ggcorrplot)
library(RColorBrewer)

setwd('/R_plot/fig2/c')
mydata = read.csv('organ_corela1.csv',row.names=1)

orRd_colors <- brewer.pal(9, "OrRd") # set the color style
ggcorrplot(mydata, hc.order = FALSE, type = "upper",
           outline.col = "gray")+ theme(panel.grid = element_blank())+
  scale_fill_gradientn(colors = orRd_colors)


#fig2 d SEM
library(lavaan)
library(semTools)
library(semPlot)

setwd('/R_plot/fig2/d')
basic_data<-read.csv("overall_inner1.csv") 

basic_model<-'  #set the model based on Tetrad
WM~GM
pancreas~liver+kidney+bonesub
kidney~liver+bonesub
heart~GM+WM+bonesub+liver
bonesub~GM
'

cfa_fit<-cfa(basic_model,basic_data) # sonfirmatory factor analysis
summary(cfa_fit,standardized=TRUE) # displays the overall summary


#fig3 wordcloud
library("wordcloud")
library("RColorBrewer")
library("ggplot2")

data<-read.csv("/R_plot/fig3/pancreas.csv")

p <- wordcloud(words = data$name, 
               freq = data$freq,
               min.freq = -15, 
               scale=c(1,0.05),
               max.words = 200,
               random.order = FALSE, 
               rot.per = 0.35)


#fig4 a
library("gplots")
library(pheatmap)
library(reshape2)

setwd('/R_plot/fig4/Survival')
mydata = read.csv('HR.csv',row.names=1)
data_matrix <- as.matrix(mydata)

p = read.csv('P.csv',row.names=1)
p_matrix <- as.matrix(p)
significance <- ifelse(p < 0.05, "*", "") 

# Generate a unique breakpoint
bk <- unique(c(seq(0.5, 1, by = 0.01), seq(1, 2, by = 0.01)))

# Create a color map
color_palette <- c(colorRampPalette(colors = c("#4E7DB7", "white"))(length(bk) * 2.2 / 7), 
                   colorRampPalette(colors = c("white", "#D63027"))(length(bk) * 4.8 / 7))

# Make sure that the length of color_palette matches the length of bk
color_palette <- color_palette[1:length(bk)]

pp <- pheatmap(mydata,
               color = color_palette,
               breaks = bk,
               angle_col = 90,
               cluster_cols = F,
               cluster_rows = F,
               border_color = "gray", 
               border_width = 0.001,
               display_numbers = significance # label "*"
)


#fig4 b
library(grid)
library(forestploter)
library(flextable)

setwd('/R_plot/fig4/single_continue')
dt <- read.csv('single.csv',,sep=',',header=TRUE)

dt$se <- 0.5
dt$` ` <- paste(rep(" ", 20), collapse = " ") #30 for multiple

dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$HR, dt$LOW, dt$UP))
p <- forest(dt[,c(1:2, 12:13)],
            est = dt$HR,       
            lower = dt$LOW, # lower limit of confidence interval
            upper = dt$UP, # upper confidence interval
            sizes = dt$se,     
            ci_column = 3,   # Select the column to draw the forest plot
            ref_line = 1,
            xlim = c(0.25, 2.4),
            ticks_at = c(0.5, 1, 2))


#fig4 c
library(grid)
library(forestploter)
library(flextable)

setwd('/R_plot/fig4/multi')
dt <- read.csv('multi.csv',,sep=',',header=TRUE)

dt$se <- 0.5
dt$` ` <- paste(rep(" ", 30), collapse = " ") #20

dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$HR, dt$LOW, dt$UP))

p <- forest(dt[,c(1, 11:12)],
            est = dt$HR,      
            lower = dt$LOW, # lower limit of confidence interval
            upper = dt$UP,  # upper confidence interval
            sizes = dt$se,    
            ci_column = 2,  # Select the column to draw the forest plot
            ref_line = 1,
            xlim = c(0.9, 2),
            ticks_at = c(1, 1.5,2))


#fig5 a
library(ggplot2)
library(ggrepel)
#Part1 Draw the scatter plot
#upload the data
setwd('/R_plot/fig5/MetabolicAsso')
data = read.csv('data4_Xhigh3.csv')
mydata = read.csv('mydata.csv')
pwas_blood_meta_sigtop10 = read.csv('pwas_blood_meta_sigtop10_adjp_2.csv')

y_max <- max(data$Log10P, na.rm = TRUE) 
y_min <- min(data$Log10P, na.rm = TRUE) -15

top_margin=0
right_margin=0
bottom_margin=0
left_margin=0

p = ggplot(data, aes(x=id, y=Log10P)) + coord_cartesian(xlim = c(0,2200),ylim = c(y_min, y_max))+
  geom_point(aes(size = Beta_abs, fill = factor(color)), color='black', alpha = 0.9, shape=21,stroke=NA) +
  scale_size_continuous(range = c(1, 4))+
  scale_fill_manual(values=c("Amino acids" = "#1F77B4", "Apolipoproteins" = "#AEC7E8", "Not significant" = '#C7C7C7',"Bone and joint" = "#FF7F0E","Cholesterol" = "#FFBB78",
                             "Endocrine" = "#2CA02C","Fatty acids" = "#98DF8A","Fluid balance" = "#FF9896","Glycolysis related metabolites" = "#9467BD","Immunometabolic" = "#C5B0D5",
                             "Inflammation" = "#8C564B","Ketone bodies" = "#C49C94","Lipoprotein and lipid" = "#E64B35","Liver function" = "#E377C2","Platelet" = "#FFD94A",
                             "Red blood cell" = "#DBDB8D","Renal function" = "#17BECF","White blood cell" = "#9EDAE5"))+  # set the color
  guides(fill = guide_legend(override.aes = list(size = 5)))+
  scale_x_continuous(expand=c(0, 0),label =mydata$groups, breaks= mydata$label) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(floor(y_min / 50) * 50, ceiling(y_max / 50) * 50, by = 50)) +
  theme_bw() +
  theme(text=element_text(),   # set the theme style
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.line = element_line(linewidth=1, colour = "black"),
        axis.text.x = element_text(color="black", size= 13,angle=45,hjust=1),
        axis.text.y = element_text(color="black", size=15),
        axis.title.y=element_text(color="black", size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(top_margin,right_margin,bottom_margin,left_margin),'lines')
  ) +
  geom_rect(aes(xmin = mydata$xstart[1],xmax = mydata$xend[1] , ymin = mydata$ystart[1] , ymax = mydata$yend[1]) ,fill="#969696",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[2],xmax = mydata$xend[2] , ymin = mydata$ystart[2] , ymax = mydata$yend[2]) ,fill="#4CAD4E",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[3],xmax = mydata$xend[3] , ymin = mydata$ystart[3] , ymax = mydata$yend[3]) ,fill="#AF72A3",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[4],xmax = mydata$xend[4] , ymin = mydata$ystart[4] , ymax = mydata$yend[4]) ,fill="#E88D91",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[5],xmax = mydata$xend[5] , ymin = mydata$ystart[5] , ymax = mydata$yend[5]) ,fill="#4189C4",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[6],xmax = mydata$xend[6] , ymin = mydata$ystart[6] , ymax = mydata$yend[6]) ,fill="#B45100",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[7],xmax = mydata$xend[7] , ymin = mydata$ystart[7] , ymax = mydata$yend[7]) ,fill="#E0BB90",alpha = 0.6)+
  
  geom_label_repel(data=pwas_blood_meta_sigtop10,aes(label=description), size= 4,
                   segment.color = "#000000", segment.size = 0.2,segment.alpha=1, force=10, box.padding =  unit(0.7, "lines"), hjust = 0.3,max.overlaps = Inf)+
  labs(x="",y="-Log10(p-value)", size="abs(Beta value)", fill= "Beta value")  + 
  theme(legend.position = "top")  # set the position of legend

#Part2 Draw the pie chart
library(ggplot2)
library(ggrepel)
library(dplyr)
#upload the data
setwd('/R_plot/fig5/MetabolicAsso/piePlot')
data = read.csv('Pancreas.csv')
#set the color
custom_color <- c("Amino acids" = "#925e9f", "Apolipoproteins" = "#BDB5E1", "nonsignificant" = '#dadada',"Bone and joint" = "#CC625F","Cholesterol" = "#F07874",
                  "Endocrine" = "#FD9BA0","Fatty acids" = "#FEC1E2","Fluid balance" = "#FCC3B4","Glycolysis related metabolites" = "#FCECCA","Immunometabolic" = "#FFFB79",
                  "Inflammation" = "#EFD964","Ketone bodies" = "#BE9E33","Lipoprotein and lipid" = "#F49513","Liver function" = "#FFD377","Platelet" = "#C0E56F",
                  "Red blood cell" = "#A0FCE8","Renal function" = "#92DDF8","White blood cell" = "#8EB4D0") # 实心八边形
#draw the plot
p <- ggplot(data, aes(x = "", y = organ, fill = color)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  scale_fill_manual(values = custom_color) +
  theme_void()


#fig5 b
library("ggplot2")
library("dplyr")
library("ggrepel")

setwd('/R_plot/fig5/ProteomicAsso')
plotdata = read.csv('imputated_median_scaled_bl_intersection_Pancreas_allremoved_result_sigResultBeta.csv')

colnames(plotdata) <- c('name','beta','p','adjp') # choose the column to draw

plotdata <- plotdata %>%
  arrange(name, desc(beta)) %>%
  distinct(name,.keep_all = TRUE)
plotdata <- as_tibble(plotdata)
plotdata$log10FDR <- -log10(plotdata$p)
plotdata$group <- case_when(plotdata$beta > 0 & plotdata$adjp < 0.05 ~ "pos", # set positive points
                            plotdata$beta < 0 & plotdata$adjp < 0.05 ~ "neg", # set negative points
                            abs(plotdata$beta) == 0 ~ "None",
                            plotdata$adjp >= 0.05 ~ "None")
top10sig<- filter(plotdata,group!= "None") %>% distinct(name,.keep_all = T) %>% top_n(10,abs(log10FDR))
plotdata$group <- factor(plotdata$group,levels= c( "pos", "neg", "None"),ordered= T)

p = ggplot(data=plotdata,aes(beta,log10FDR,color=group))+
  geom_point(size=2.5,stroke=NA)+
  scale_colour_manual(name= "",values=alpha(c( "#CB4301", "#1065C0", "gray80"),0.9))+
  geom_text_repel(data=top10sig, aes(beta, log10FDR, label=name), max.overlaps=20)+
  theme_classic()+
  geom_hline(yintercept = c(-log10(0.05)), linewidth= 0.7,color= "grey",lty= "dashed")+
  geom_vline(xintercept = c(0),size= 0.7, color= "grey", lty= "dashed")+
  labs(y = "-log10P")


#fig5 c
library(ggpubr)

setwd('/R_plot/fig5/ProteomicAsso/c_lillipop')
data = read.csv('tissue.csv')

#Custom sort vector
custom_order <- c('GM_Brain',
                  'GM_Esophagus',
                  'GM_Nerve',
                  'GM_Prostate',
                  'GM_Pituitary',
                  'WM_Brain',
                  'WM_Esophagus',
                  'WM_Stomach',
                  'WM_Colon',
                  'WM_Salivary_Gland',
                  'heart_Breast',
                  'heart_Kidney',
                  'heart_Adipose_Tissue',
                  'heart_Lung',
                  'heart_Testis',
                  'bonesub_Heart',
                  'bonesub_Muscle',
                  'bonesub_Esophagus',
                  'bonesub_Nerve',
                  'bonesub_Colon',
                  'kidney_Kidney',
                  'kidney_Pancreas',
                  'kidney_Testis',
                  'kidney_Brain',
                  'kidney_Pituitary',
                  'liver_Liver',
                  'liver_Kidney',
                  'liver_Adipose_Tissue',
                  'liver_Breast',
                  'liver_Stomach',
                  'pancreas_Pancreas',
                  'pancreas_Stomach',
                  'pancreas_Liver',
                  'pancreas_Kidney',
                  'pancreas_Ovary') 
# Convert the "tissue" columns to factors and set a custom order
data$tissue <- factor(data$tissue, levels = custom_order)

p <- ggplot(data, aes(x=tissue, y=value, colour = group)) +
  geom_segment(aes(x=tissue, xend=tissue, y=0, yend=value)) +
  geom_point(size=4,aes(colour = group)) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#706d94","#ecd09c","#a3c8a4","#b96570"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#706d94","#ecd09c","#a3c8a4","#b96570"))+
  theme_light()+
  theme(
    panel.grid.major.y = element_blank(),   # Remove the main horizontal gridlines
    panel.grid.minor.y = element_blank(),   # Remove secondary horizontal gridlines
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(colour = "black") # Add axis
  )
p + geom_hline(yintercept = -log10(0.00173), linetype = "dashed", color = "black")

# fig 6
library(ggplot2)
library(ggrepel)

#upload the data
setwd('/R_plot/fig6')
data = read.csv('data_Xhigh3.csv')
mydata = read.csv('mydata.csv')
pwas_blood_meta_sigtop10 = read.csv('pwas_blood_meta_sigtop10_adjp.csv')

y_max <- max(data$Log10P, na.rm = TRUE) 
y_min <- min(data$Log10P, na.rm = TRUE) -15

top_margin=0
right_margin=0
bottom_margin=0
left_margin=0

#plot the figure , 
p = ggplot(data, aes(x=id, y=Log10P)) + coord_cartesian(xlim = c(0,1200),ylim = c(y_min, y_max))+
  geom_point(aes(size = Beta_abs, fill = factor(color)), color='black', alpha = 0.9, shape=21,stroke=NA) +
  scale_fill_manual(values=c("Cognitive function" = "#F57C6E","Physical measures" = "#F2B56F", "nonsignificant" = '#dadada',"Health and medical history" = "#84C3B7",
                             "Lifestyle and environment" = "#71B7ED","Psychosocial factors" = "#FAE69E","Sociodemographics" = "#B8AEEB",
                             "Early life and family history factors" = "#F2A7DA"))+   # set the color
  guides(fill = guide_legend(override.aes = list(size = 5)))+
  scale_x_continuous(expand=c(0, 0),label =mydata$groups, breaks= mydata$label) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(floor(y_min / 50) * 50, ceiling(y_max / 50) * 50, by = 50)) +
  theme_bw() +
  theme(text=element_text(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.line = element_line(linewidth=1, colour = "black"),
        axis.text.x = element_text(color="black", size= 13,angle=45,hjust=1),
        axis.text.y = element_text(color="black", size=15),
        axis.title.y=element_text(color="black", size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(top_margin,right_margin,bottom_margin,left_margin),'lines')
  ) +
  geom_rect(aes(xmin = mydata$xstart[1],xmax = mydata$xend[1] , ymin = mydata$ystart[1] , ymax = mydata$yend[1]) ,fill="#969696",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[2],xmax = mydata$xend[2] , ymin = mydata$ystart[2] , ymax = mydata$yend[2]) ,fill="#4CAD4E",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[3],xmax = mydata$xend[3] , ymin = mydata$ystart[3] , ymax = mydata$yend[3]) ,fill="#AF72A3",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[4],xmax = mydata$xend[4] , ymin = mydata$ystart[4] , ymax = mydata$yend[4]) ,fill="#E88D91",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[5],xmax = mydata$xend[5] , ymin = mydata$ystart[5] , ymax = mydata$yend[5]) ,fill="#4189C4",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[6],xmax = mydata$xend[6] , ymin = mydata$ystart[6] , ymax = mydata$yend[6]) ,fill="#B45100",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[7],xmax = mydata$xend[7] , ymin = mydata$ystart[7] , ymax = mydata$yend[7]) ,fill="#E0BB90",alpha = 0.6)+
  geom_rect(aes(xmin = mydata$xstart[8],xmax = mydata$xend[8] , ymin = mydata$ystart[8] , ymax = mydata$yend[8]) ,fill="#7DDFD7",alpha = 0.6)+
  geom_label_repel(data=pwas_blood_meta_sigtop10,aes(label=description), size= 4,
                   segment.color = "#000000", segment.size = 0.2,segment.alpha=1, force=10, box.padding =  unit(0.7, "lines"), hjust = 0.3,max.overlaps = Inf)+
  labs(x="",y="-Log10(p-value)", size="abs(Beta value)", fill= "Beta value")  + 
  theme(legend.position = "top") # set the position of legend


# fig 7 a
library("gplots")
library(pheatmap)
library(reshape2)

setwd('/R_plot/fig7')
mydata = read.csv('HR.csv',row.names=1)

mydata1 <- as.data.frame(apply(mydata, 2, function(x) ifelse(x < 0.5, 0, x)))
bk <- c(seq(0.5,0.69,by=0.01),seq(0.7,0.85,by=0.01))
p <- pheatmap(mydata1,
              #color =colorRampPalette(colors=c("white","#D63027"))(100),
              #color = c(colorRampPalette(colors = c("white","#D63027"))),
              color = c(colorRampPalette(colors = c("white","#ffe8b8"))(length(bk)*1/5),colorRampPalette(colors = c("#ffe8b8","#D63027"))(length(bk)*4/5)),
              angle_col = 90,
              breaks = bk,
              # scale="row",
              # border="white", 
              cluster_cols = F,
              cluster_rows = F
)

# fig 7 b
library(glmnet)
library(caret)
library(ggplot2)
library(gridExtra)
library(pROC)
library(pracma)

setwd('/R_plot/fig7')
mydata = read.csv('E2.csv')
mydata1 = read.csv('F1.csv')
mydata2 = read.csv('I4.csv')
mydata3 = read.csv('X2.csv')
mydata4 = read.csv('E1.csv')

auc <- trapz(mydata$fpr, mydata$tpr)
input = mydata
auc1 <- trapz(mydata1$fpr, mydata1$tpr)
input1 = mydata1
auc2 <- trapz(mydata2$fpr, mydata2$tpr)
input2 = mydata2
auc3 <- trapz(mydata3$fpr, mydata3$tpr)
input3 = mydata3
auc4 <- trapz(mydata4$fpr, mydata4$tpr)
input4 = mydata4
source('ztheme_ROC.R')

color <- c("Obesity" = '#F86DF0',
           "Dementia" = '#9FDB9F', 
           "Heart failure" = '#797CE1', 
           "Nervous system" = '#acc9e9',
           "Diabetes" = '#E68b81')

p <- ggplot() +
  geom_line(data = input, aes(x = fpr, y = tpr, color = "Obesity"), size = 1) +
  geom_line(data = input1, aes(x = fpr, y = tpr, color = "Dementia"), size = 1) +
  geom_line(data = input2, aes(x = fpr, y = tpr, color = "Heart failure"), size = 1) +
  geom_line(data = input3, aes(x = fpr, y = tpr, color = "Nervous system"), size = 1) +
  geom_line(data = input4, aes(x = fpr, y = tpr, color = "Diabetes"), size = 1) +
  geom_abline(intercept = 0, size = 1, alpha = 0.5, color = 'gray', slope = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_color_manual(values = color) +
  labs(x = "False Positive Rate", 
       y = "True Positive Rate", 
       color = "Legend") + 
  z_theme() 

ggsave("ALL.pdf", plot = p, width = 5.2, height = 4) # save the plot