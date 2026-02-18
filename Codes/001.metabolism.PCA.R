rm(list = ls())
metainfo = openxlsx::read.xlsx("../Data/000.clinical.process.addCell.xlsx")

dat = read.csv("./org_mix/all_metab_data.txt",sep = "\t")

dat_meta = dat[ , 23:ncol(dat )]
row.names(dat_meta) <- dat$metab_id
dat_meta = dat_meta[ , as.character( metainfo$uID ) ]


dat_meta = apply(dat_meta,2,function(x){10^6 *x/sum(x) } )


dat_meta = dat_meta[ apply( dat_meta,1,function(x){sum(x>0)  > 5  } ) ,]
dat_rela_tmp = dat_meta

dat_rela_tmp = apply(dat_rela_tmp, 1, function(x){ sd(x)/mean(x)   }   )
dat_rela_tmp = sort(dat_rela_tmp , decreasing = T)
dat_meta = dat_meta[   names(dat_rela_tmp) %>% head( length(dat_rela_tmp) *0.75  ), ]


select_time = c("001_BL","002_D1","003_D14","004_M1","005_M2","006_M3","009_M6","013_M12")
metainfo = metainfo[metainfo$timeIDcontiue %in% select_time ,]
dat_meta = dat_meta[  ,metainfo$uID ]
dat_meta = dat_meta[ apply( dat_meta,1,function(x){sum(x>0)  > 5  } ) ,]


## —— 可选：对强度做 log 转换（推荐）——
## 代谢组强度通常右偏，先 log1p 有助于稳定方差
dat_meta_log <- log1p(dat_meta)

## —— 三种缩放方法：CTR / UV / PAR ——
scale_ctr <- function(X){               # 仅中心化
  scale(X, center = TRUE, scale = FALSE)
}
scale_uv  <- function(X){               # autoscaling：/SD
  scale(X, center = TRUE, scale = TRUE)
}
scale_par <- function(X){               # Pareto：/sqrt(SD)
  sds <- apply(X, 2, sd, na.rm = TRUE)
  scale(X, center = TRUE, scale = sqrt(sds))
}

## 为了与 prcomp 一致：行=样本，列=代谢物
X <- t(dat_meta_log)   # 你前面已经做过 t()，这里对 log 后的矩阵再做一次

## 跑 PCA 并绘图的小函数
library(factoextra)
plot_pca <- function(Xscaled, title_txt){
  pcs <- prcomp(Xscaled, center = FALSE, scale. = FALSE)  # 已手动缩放过，这里都设 FALSE
  p <- fviz_pca_ind(
    pcs,palette = c("grey20","#feb24c","red", "#7FC97F", "#BEAED4", "#FDC086", "#c51b8a", "#386CB0" ),
    geom = "point", ellipse.level = 0.95,
    # pointshape = c(1),  
    col.ind = metainfo$timeIDcontiue, addEllipses = TRUE,
    ellipse.type = "confidence", legend.title = "Groups", repel = TRUE
  ) +
    ggplot2::theme_bw(12) +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text = ggplot2::element_text(family = "ArialMT"),
      axis.title = ggplot2::element_text(family = "ArialMT"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey87")
    ) + scale_shape_manual( values = c(1,2,3,4,5,6,7,8) ) + scale_size_manual( values = c(10,1) ) +
    ggplot2::ggtitle(title_txt)
  return(p)
}

## 生成三种缩放的 PCA 图
p_log <- plot_pca(t(dat_meta_log), "PCA-log")
p_ctr <- plot_pca(scale_ctr(X), "PCA-CTR")
p_uv  <- plot_pca(scale_uv(X),  "PCA-UV")
p_par <- plot_pca(scale_par(X), "PCA-PAR")

ggsave("./metabolism/001.uv_pva.pdf",p_uv,width = 5.79,height = 4.21)
