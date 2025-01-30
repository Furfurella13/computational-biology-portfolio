#Данные: нормализованная экспрессия генов при мелкоклеточном раке легкого (DOI: 10.1038/nature14664).
# Определить дифференциальные гены между пациентами IV и I стадии. Метод – любой, применимый к имеющимся данным. Критерии определения гена как дифференциального – p-value adjusted<0.05, |log2FoldChange|>2. Визуализировать в формате volcano plot( предпочтительно через ggplot2).
# Для каждого семейства генов посчитать среднюю экспрессию(ср.геометрическое на пациента), определить корреляцию со стадией  и p-value для корреляции. Выполнить поправку на множественную проверку гипотез, метод - fdr. Сохранить в файл .xlsx (формат - семейство генов, средние у пациентов, коэффициент корреляции, p-value, p-adjusted). Построить хитмап со средними, разметить пациентов согласно стадиям с помощью side-bar(предпочтительно pheatmap).  При наличии семейств генов, значимо коррелирующих со стадией процесса(p-adjusted<0.05), для топ-5 семейств с отрицательной корреляцией и топ-5 с положительной корреляцией, визуализировать распределение в виде боксплотов для каждой стадии(предпочтителен boxplot  - базовая графика). 
# Все рисунки надо сохранить в формате tiff, шириной в 7,5 дюйма, с разрешением >=300 dpi. Шрифты должны быть читаемы при размещении на странице А4.

# Загрузка библиотек
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("edgeR")
BiocManager::install ("limma")

#library(DESeq2)
library(tidyverse)
library(pheatmap)
#library(edgeR)
library(ggplot2)
library(tibble)
library(tidyr)
setwd("C:\\Users\\nasty\\Desktop\\Магистратура\\R\\Для зачета")

file_path <- "data_clinical_patient.txt"
#file_path <- "C:/Users/White knight/Downloads/Telegram Desktop/data_clinical_patient.txt"
data <- read.delim(file_path, header = TRUE, sep = "\t", dec = ".")
#head(data, n=5)
file_path2 <- "data_mrna_seq_rpkm.txt"
#file_path2 <- "C:/Users/White knight/Downloads/Telegram Desktop/data_mrna_seq_rpkm.txt"
data2 <- read.delim(file_path2, header = TRUE, sep = "\t", dec = ".")
data3 <- data2 %>% select(!("Entrez_Gene_Id"))
#data2 <- data2 %>% group_by(Hugo_Symbol) %>% summarise_at(vars(-group_cols()), mean, na.rm = TRUE)
# Удаление дубликатов и вычисление среднего значения данных по каждому гену
data3 <- aggregate(. ~ Hugo_Symbol, data = data3, FUN = mean, na.rm = TRUE)
#head(data2, n=5)
#data2 <- data2 [-which(duplicated(data2$Hugo_Symbol)),]
rownames(data3) <- data3$Hugo_Symbol
data3 <- data3 %>% select(!("Hugo_Symbol"))

head(data3, n = 3, c(3))
# # Перевести широкий формат в длинный
# data_long <- gather(data2, key = "PATIENT_ID", value = "value", -c(1:17))
# 
# # Просмотреть результат
# head(data_long)
# 
# df <- data.frame(data_long$Hugo_Symbol)
# df$PATIENT_ID <- data_long$PATIENT_ID
# df$v <- rep(1, nrow(df))
# 
# dfq <- df %>%  pivot_wider( names_from = PATIENT_ID, values_from=v, values_fn = sum, values_fill = 0)

# Посмотрим на уникальные значения в столбце UICC_TUMOR_STAGE
unique_values <- unique(data$UICC.Tumor.Stage)
print(unique_values)

# Удаление букв "a" и "b" из столбца "UICC.Tumor.Stage"
data <- mutate(data, UICC.Tumor.Stage = gsub("[abAB]", "", UICC.Tumor.Stage)) 

# Выберем PATIENT_ID для стадий I и IV
selected_patient_ids <- data$X.Patient.Identifier[data$UICC.Tumor.Stage %in% c("I", "IV")]

# Выведем выбранные PATIENT_ID
print(selected_patient_ids)

#data2 <- data2 %>% select(!("Entrez_Gene_Id"))

need_patients <- intersect(selected_patient_ids, colnames(data3))

# Выберем нужные PATIENT_ID
selected_patients <- data3[,need_patients] 
#head(selected_patients, n = 5)
# Размеры selected_patients
print(dim(selected_patients))

stage_vector <- data$UICC.Tumor.Stage [match(need_patients,data$X.Patient.Identifier)]
print(stage_vector)

library(limma)

design <- cbind(Grp1 = 1, Grp2vs1 = as.numeric(stage_vector=="IV"))

# Ordinary fit
fit <- lmFit(selected_patients, design)
fit <- eBayes(fit)
dim(fit)
res <- topTable(fit, coef = 2,number = nrow(fit)) 

#head(res, n=5)
#head(res)
#dim(res)
#unique(stage_vector)
#table(stage_vector)


# Определение дифференциальных генов
#diff_genes <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 2, ]
diff_genes <- res[res$P.Value < 0.05 & abs(res$logFC) > 2, ]
#print(sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 2))

# Визуализация в вулкане
volcano_plot <- ggplot(diff_genes, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(abs(logFC) > 2, "red", "black")), alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray") +
  labs(title = "Volcano Plot of Differential Expression",
       x = "log2 Fold Change",
       y = "-log10(P.Value)")

# Показать вулкан
print(volcano_plot)

data_1 <- read.table('C:/Users/nasty/Desktop/Магистратура/R/Для зачета/data_mrna_seq_rpkm.txt', sep = '\t')
raw_counts <- as.data.frame(lapply(data_1[-1, -c(1, 2)], as.numeric))
raw_counts <- raw_counts[-which(duplicated(data_1[-1, 1])), ]
colnames(raw_counts) <- data_1[1, -c(1, 2)]
rownames(raw_counts) <- data_1[-1, 1][-which(duplicated(data_1[-1, 1]))]

data_2 <- read.table('C:/Users/nasty/Desktop/Магистратура/R/Для зачета/data_clinical_patient.txt', sep = '\t')
meta_data <- data_2[-1, c(1, 7)]
meta_data <- meta_data[which(data_2[-1, 1] %in% data_1[1, -c(1, 2)]),]

stage <- c()
for (i in 1:length(meta_data[, 2])){
  if (i %in% grep('IV', meta_data[, 2])){
    stage <- c(stage, '4')
  }
  else if (i %in% grep('III', meta_data[, 2])){
    stage <- c(stage, '3')
  } 
  else if (i %in% grep('II', meta_data[, 2])){
    stage <- c(stage, '2')
  }
  else if (i %in% grep('I', meta_data[, 2])){
    stage <- c(stage, '1')
  }
  else {
    stage <- c(stage, NaN)
  }
}
meta_data <- as.data.frame(cbind(meta_data[, -2], as.numeric(stage))) 
colnames(meta_data) <- c('ID', 'STAGE')

stage_1 <- meta_data[which(meta_data$STAGE == 1), 1]
stage_2 <- meta_data[which(meta_data$STAGE == 2), 1]
stage_3 <- meta_data[which(meta_data$STAGE == 3), 1]
stage_4 <- meta_data[which(meta_data$STAGE == 4), 1]

data_stage_1 <- raw_counts[, which(colnames(raw_counts) %in% stage_1)]
data_stage_2 <- raw_counts[, which(colnames(raw_counts) %in% stage_2)]
data_stage_3 <- raw_counts[, which(colnames(raw_counts) %in% stage_3)]
data_stage_4 <- raw_counts[, which(colnames(raw_counts) %in% stage_4)]

data_3 <- read.table('C:/Users/nasty/Desktop/Магистратура/R/Для зачета/hgnc_complete_set_2017.txt', sep = '\t',  fill=TRUE)
data_family <- data_3[-1, c(2, 13)]
data_family <- data_family[-which(nchar(data_family$V13) == 0),]

nom_family <- which(rownames(raw_counts) %in% data_family$V2)
gene_family <- c()
for (i in 1:length(nom_family)){
  gene_family[i] <- data_family[which(data_family$V2 %in% rownames(raw_counts[nom_family, ])[i]), 2]
}

fam_st_1 <- cbind(data_stage_1[nom_family, ], as.factor(gene_family))
fam_st_2 <- cbind(data_stage_2[nom_family, ], as.factor(gene_family))
fam_st_3 <- cbind(data_stage_3[nom_family, ], as.factor(gene_family))
fam_st_4 <- cbind(data_stage_4[nom_family, ], as.factor(gene_family))

log_mean <- function(x){
  exp(mean(log(x[x > 0])))
}

f1 <- fam_st_1 %>%
  group_by(fam_st_1$`as.factor(gene_family)`) %>%
  summarise_if(is.numeric, log_mean)

f2 <- fam_st_2 %>%
  group_by(fam_st_2$`as.factor(gene_family)`) %>%
  summarise_if(is.numeric, log_mean)

f3 <- fam_st_3 %>%
  group_by(fam_st_3$`as.factor(gene_family)`) %>%
  summarise_if(is.numeric, log_mean)

f4 <- fam_st_4 %>%
  group_by(fam_st_4$`as.factor(gene_family)`) %>%
  summarise_if(is.numeric, log_mean)

p_values_family <- c()
coef_cor <- c()
for (i in 1:nrow(f1)){
  y <- c(rep(1, length(f1[i, -1])), rep(2, length(f2[i, -1])),
         rep(3, length(f3[i, -1])), rep(4, length(f4[i, -1]))) 
  x <- c(as.numeric(f1[i, -1]), as.numeric(f2[i, -1]),
         as.numeric(f3[i, -1]), as.numeric(f4[i, -1]))
  if (length(x) - sum(is.na(x)) > 1){
    cor_sp <- cor.test(x, y, method = 'spearman')
    p_values_family[i] <- cor_sp$p.value
    coef_cor[i] <- as.numeric(cor_sp$estimate) 
  } else {
    p_values_family[i] <- NaN
    coef_cor[i] <- NaN
  }
}

p_adjust_family <- p.adjust(p = p_values_family, method = 'fdr')

my_result <- cbind(f1[, 1], coef_cor, p_values_family, p_adjust_family)

library(openxlsx)
write.xlsx(my_result, "my_result.xlsx",  rowNames = TRUE)

combined_data <- cbind(f1[, -1], f2[, -1], f3[, -1], f4[, -1])

md_mat <- as.matrix(combined_data) 
heat_map_data <- log(md_mat + 1) 
# Замена NA на среднее значение
heat_map_data[is.na(heat_map_data)] <- mean(heat_map_data, na.rm = TRUE)

my_col <- c()
for (i in 1:length(colnames(heat_map_data))){
  if (colnames(heat_map_data)[i] %in% colnames(data_stage_1)){
    my_col <- c(my_col, '#9ACD32')
  }
  else if (colnames(heat_map_data)[i] %in% colnames(data_stage_2)){
    my_col <- c(my_col, '#CD8500')
  }
  else if (colnames(heat_map_data)[i] %in% colnames(data_stage_3)){
    my_col <- c(my_col, '#CD3700')
  }
  else if (colnames(heat_map_data)[i] %in% colnames(data_stage_4)){
    my_col <- c(my_col, '#CD0000')
  } 
  else{
    my_col <- c(my_col, 'black')
  }
}  

heatmap(heat_map_data, ColSideColors = my_col, Colv = NA)

# Выбор топ-5 семейств с отрицательной корреляцией
top_negative <- my_result[my_result$coef_cor < 0, ]
top_negative <- top_negative[order(top_negative$p_adjust_family), ][1:5, ]

# Выбор топ-5 семейств с положительной корреляцией
top_positive <- my_result[my_result$coef_cor > 0, ]
top_positive <- top_positive[order(top_positive$p_adjust_family), ][1:5, ]

# Объединение топ-5 семейств в один список
top_gene_families <- rbind(top_negative, top_positive)

# Построение боксплотов для каждой стадии рака


library(ragg)
# Установка параметров графического устройства
par(mar = c(5, 4, 4, 5), cex.main = 1.2, cex.lab = 1.0, cex.axis = 0.8)

# Создание файла TIFF
agg_tiff(filename = "Volcano_Plot.tiff", width = 7.5, height = 7.5, units = "in", res = 300)
plot(volcano_plot, main = "Volcano Plot")
dev.off()

# Установка параметров графического устройства
par(mar = c(5, 4, 4, 5))

# Создание файла TIFF
agg_tiff(filename = "Heat_map_stage.tiff", width = 7.5, height = 7.5, 
         units="in", res=300)
heatmap(heat_map_data, ColSideColors = my_col, Colv = NA)
dev.off()