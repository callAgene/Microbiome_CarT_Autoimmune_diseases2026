# library(edgeR)
rm(list = ls())

metainfo = openxlsx::read.xlsx("./Data/000.clinical.process.addCell.xlsx")
metainfo = metainfo[  metainfo$timeIDcontiue %in% c("001_BL","002_D1","003_D14","004_M1","005_M2","006_M3","009_M6","0013_M12")  , ]

# 读入数据
dat <- read.delim( gzfile( "./Data/003.metagenome/001.KO_unstratified.tsv.gz"), row.names=1, check.names=FALSE)

# 不移除非功能行（如 UNMAPPED、UNINTEGRATED）
colnames(dat ) <- sub("_genefamilies","",colnames( dat ))

dat <- dat[  ! grepl( "g__",row.names(dat) ) ,]
dat = dat[  , metainfo$metaID ]

row.names(metainfo) <- metainfo$metaID

dat = apply( dat,2,function(x) { 10^6*x/sum(x) }  )
dat = dat[ apply(dat,1,function(x) { sum(x>0) > length(x) *0.2  } ) ,]

pseudo <- 1e-6
library(limma)

y <- log2(dat + 1e-6)    


metainfo_select <- metainfo



# 对齐样本顺序
stopifnot(all(colnames(y) == rownames(metainfo_select)))

# 2) 分组因子（按你的时间顺序设 levels）
group <- factor(
  metainfo_select$timeIDcontiue,
  levels = c("001_BL","002_D1","003_D14","004_M1","005_M2","006_M3","009_M6")
)
levels(group) <- make.names(levels(group))

# 3) 设计矩阵（无截距）
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 5) 拟合 + 趋势 eBayes（适合 log 强度数据）
fit <- lmFit( y, design)                 
fit <- eBayes(fit, trend = TRUE)

# 6) 定义组间对比（示例：D1 vs BL、D14 vs BL、M1 vs BL）
cont <- makeContrasts(
  D1_vs_BL  = `X002_D1` - `X001_BL`,
  D14_vs_BL = `X003_D14` - `X001_BL`,
  M1_vs_BL  = `X004_M1` - `X001_BL`,
  M2_vs_BL  = `X005_M2` - `X001_BL`,
  M3_vs_BL  = `X006_M3` - `X001_BL`,
  M6_vs_BL  = `X009_M6` - `X001_BL`,
  M1_vs_D14  = `X004_M1` - `X003_D14`,
  M12_vs_BL  = `X004_M1` - `X003_D14`,
  levels = design
)


fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE)


########################################################### 
# 7) 取某个对比的结果（例如 D14 vs BL）
result=list()

input_com = "D1_vs_BL"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt


input_com = "D14_vs_BL"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt


input_com = "M1_vs_BL"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt


input_com = "M2_vs_BL"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt


input_com = "M3_vs_BL"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt

input_com = "M6_vs_BL"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt


input_com = "M1_vs_D14"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt


input_com = "M12_vs_BL"
tt <- topTable(fit2, coef = input_com, n = Inf, adjust.method = "BH")
tt$pathway <- rownames(tt)   # 给结果加上特征名
tt$type=input_com
result[[input_com]] = tt

result[["cpm_dat"]] =  data.frame(pathway = row.names(dat),dat)
result[["metainfo"]] = metainfo



saveRDS(result,"./002.metagenome/001.KO.pathway.RDS")

openxlsx::write.xlsx(result,"./002.metagenome/001.DF.KO.xlsx")


kegg_dat = readRDS("./Data/KEGG/ko_long.RDS")
kegg_dat = kegg_dat[ ! kegg_dat$A %in% c("A09160 Human Diseases") ,]
kegg_dat = kegg_dat[  kegg_dat$A %in% c("A09100 Metabolism") ,]

kegg_dat_name = kegg_dat[  , c( "C_code","C_name" ) ] %>% unique(.)
colnames( kegg_dat_name  ) <- c( "term" ,"name" )

kegg_dat_select = kegg_dat[,c("C_code","K")] 
colnames(kegg_dat_select) <- c("term" ,"gene")


library(clusterProfiler)

# example for D14_vs_BL : 
select_dat = result[["D14_vs_BL"]]
select_dat_gene = select_dat$logFC
names(select_dat_gene) <-  select_dat$pathway 
select_dat_gene = sort(select_dat_gene,decreasing = T)

kegg_dat_select = kegg_dat_select[ kegg_dat_select$gene %in%   row.names(select_dat),]
result_GSEA = GSEA(geneList = select_dat_gene,TERM2GENE =   kegg_dat_select,TERM2NAME =  kegg_dat_name )

result_GSEA@result


paths <- c("ko00360", "ko00121", "ko00640")#选取你需要展示的通路ID
enrichplot::gseaplot2(result_GSEA,paths, pvalue_table = TRUE)
gseaplot2(result_GSEA,paths,color = colorspace::rainbow_hcl(4),subplots=c(1,2), pvalue_table = TRUE)

