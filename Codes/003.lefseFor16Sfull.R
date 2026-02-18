
  
calculate_timeID = "003_D14"
loadtmp = readRDS("../Data/002.16Sfull/001.species.dat.Rds");
dat_rela = loadtmp$dat_rela

############
metainfo = openxlsx::read.xlsx("../Data/000.clinical.process.xlsx",sheet = "S16")


############


dat_rela_select = dat_rela_select[ rowSums(dat_rela_select) > 0 ,]

metainfo_select =  metainfo
row.names(metainfo_select) <- metainfo_select$uID
metainfo_select = metainfo_select[ colnames(dat_rela_select  ) ,]

all_patient = metainfo_select$patID %>% unique(.) 


library(lefser)

dat_name = read.csv("../Data/002.16Sfull/ASV_Taxon_Origin_asv.full.xls",sep = "\t")
dat_name = dat_name[,c(1:8)]
dat_name = dat_name[,-c(2)]

dat_name$merger = apply(dat_name , 1,paste,sep = "",collapse = "|")
dat_name = unique(dat_name)
dat_name = dat_name[ match(  row.names(dat_rela_select) ,dat_name$species  ) , ]


row.names(dat_name) <- dat_name$species
table( row.names(dat_name) == row.names( dat_rela_select ) )

## Split data tables
counts <- as.data.frame(dat_rela_select)
coldata <- as(metainfo_select, "data.frame")

## Create a SummarizedExperiment object
objectD = SummarizedExperiment(assays = list(counts = counts),
                               colData = coldata,
                               rowData = dat_name[,-8])




objectD_select  = objectD[  ,objectD$timeIDcontiue %in% c("001_BL","002_D1") ]

set.seed(1234)
res <- lefser(objectD_select, # relative abundance only with terminal nodes
              groupCol = "timeIDcontiue")
head(res)

lefserPlot(res)


all_time =  c("002_D1"  ,"003_D14", "004_M1",  "005_M2"  ,"006_M3"  ,"007_M4"  ,"008_M5"  ,"009_M6" , "010_M8",  "011_M9",  "012_M10", "013_M12" )

result = list()
for (each in all_time ){
# each = "002_D1"

objectD_select  = objectD[  ,objectD$timeIDcontiue %in% c("001_BL",each) ]

set.seed(1234)
res <- lefser(objectD_select, # relative abundance only with terminal nodes
              groupCol = "timeIDcontiue")

if ( nrow(res)  == 0) { next }
res$group = ifelse( res$scores > 0 , each , "001_BL"  )



res$Names = factor(res$Names  , levels = as.character( res$Names  )  )
p <- ggplot(as.data.frame(res), aes( x = Names, y = scores) ) + 
        ylab("LDA SCORE (log 10)") + theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 11, face = "bold"), 
        axis.text.y = element_text(vjust = 0.7, size = 9, face = "bold"), 
        axis.text.x = element_text(vjust = 0.7, size = 9, face = "bold"), 
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold")) + 
       geom_bar(stat = "identity", aes(fill = group),width = 0.618) + scale_fill_manual(values = c("#7da1cc","#e57a77")) + 
        coord_flip()
ggsave( sprintf("./Code/001.16S/012.lefse/%s.pdf",each),p,width = 5.82,height =6.32)                 
result[[ each  ]]   =   res
}

result[["dat_name"]] = dat_name
result[["counts"]] = data.frame( row=row.names(counts) , counts )
result[["coldata"]] = coldata



openxlsx::write.xlsx( result,"./001.16S/012.lefse/total.xlsx" )
】