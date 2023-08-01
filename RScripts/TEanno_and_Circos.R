#!/usr/bin Rscript

pacman::p_load(tidyverse, circlize,data.table,kableExtra)

setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/TME117-assembly/FinalAssembly/TEAnnotation/")

TE <- fread("TME117Scaf18Hap1Chr.fasta.mod.EDTA.TEanno.gff3" , h=F)

TE <- TE[, .(V1,V4,V5,V2,V9,V3)]

colnames(TE) <- c("chr", "start", "end", "id", "info", "classification")

#View(TE[classification=="repeat_region"])

#make classification of the annotated TEs - using the classification column 
REPEATS <- TE %>% mutate(
  tmp_name = str_extract(info, "Name=[^;]*"),
  family = str_extract(tmp_name, "[^=]*$"),
  tmp_subclass = str_extract(info, "Classification=[^;]*"),
  subclass = str_extract(tmp_subclass, "[^=]*$")
)%>% 
select(chr, start, end, id, family, subclass, classification) %>%
  filter(str_detect(classification, "TIR|LTR|helitron|LINE_element|repeat_region")) %>%
  separate(subclass, c("class", "superfamily"), sep = "/") %>%
  mutate(
    order = case_when(
      str_detect(classification, "TIR") & !str_detect(class, "MITE") ~ "TIR",
      str_detect(class, "MITE") ~ "MITE",
      str_detect(classification, "helitron") ~ "Helitron",
      str_detect(classification, "LTR") ~ "LTR-RT",
      str_detect(classification, "LINE") ~ "LINE",
      str_detect(classification, "repeat_region") ~ "Unclassified repeat"
    )
)%>% mutate(
      class = case_when(
      order == "TIR" | order == "MITE" | order == "Helitron" ~ "DNA",
      order == "LTR-RT" ~ "retrotransposon",
      order == "LINE" ~ "non-LTR retrotranspson",
      order == "Unclassified repeat" ~ "Other"
    ),
    superfamily = if_else(is.na(superfamily), order, superfamily),
    superfamily = case_when(
      str_detect(superfamily, "Gypsy") ~ "Gypsy",
      str_detect(superfamily, "Copia") ~ "Copia",
      str_detect(superfamily, "Helitron") ~ "Helitron",
      str_detect(superfamily, "LINE") ~ "LINE",
      str_detect(superfamily, "unknown") ~ "unknown",
      str_detect(superfamily, "DTM|Mutator") ~ "MUDR-Mutator",
      str_detect(superfamily, "DTH|Harbinger") ~ "PIF-Harbinger",
      str_detect(superfamily, "DTA|hAT") ~ "hAT",
      str_detect(superfamily, "DTC|CACTA") ~ "CACTA",
      str_detect(superfamily, "DTT|Mariner") ~ "Tc1-Mariner"
    ),
    superfamily = if_else(is.na(superfamily), "-", superfamily)
)

#View(REPEATS)
#####################################################
##### summary stats 

# Unique families
families <- REPEATS %>%
  select(family, superfamily, order, class) %>%
  distinct(family, .keep_all = TRUE)

# Number of unique families
copy_num <- REPEATS %>%
  group_by(family) %>%
  summarise(num = n()) %>%
  arrange(desc(num)) %>%
  left_join(families, by = "family")
################################################

# Length of each annotated element
anno_start <- REPEATS %>%
  mutate(
    class = case_when(
      str_detect(order, "LTR|LINE") ~ "Class I",
      str_detect(order, "Helitron|TIR|MITE") ~ "Class II",
      str_detect(order, "Unclassified repeat") ~ "Other"
    ),
    length = end - start
  ) %>%
  select(chr,class, order, superfamily, family, id, start, end,length)

# bp annotated for each superfamily
anno_len <- anno_start %>%
  group_by(class, order, superfamily) %>%
  summarise(across(length, sum))

# Number of families in each superfamily
anno_fam <- anno_start %>%
  distinct(family, .keep_all = TRUE) %>%
  group_by(class, order, superfamily) %>%
  summarise(across(family, ~ n()))

# Summary table for annotated TEs
anno_summary <- anno_start %>%
  group_by(class, order, superfamily) %>%
  summarise(across(id, ~ n())) %>%
  left_join(anno_len, by = c("class", "order", "superfamily")) %>%
  left_join(anno_fam, by = c("class", "order", "superfamily")) %>%
  mutate(
    prop = (length / 677250621) * 100,
    prop = round(prop, digits = 2)
  ) 

write_csv(anno_summary, "TEsummaryHAp1.csv")

#remove rows with "other"
anno_sum <- anno_summary[anno_summary$order!="Unclassified repeat",]
#subset(anno_summary, order!="Unclassified repeat")
anno_sum %>%
  select(class,order,superfamily,prop) %>%
  kbl(
    col.names = c(
      "Class",
      "Family",
      "Superfamily",
      "Total %"
    ),
    align = c("l", "l", "c", "c"),
    caption = "Table: The proportion of different classes of annotated transposons  
     occupied in the genome of <i> Manihot esculenta</i> - TME 117 genome - Hap1"
  ) %>%
  kable_classic(full_width = F, html_font = "Helvetica") %>%
  kable_styling(font_size = 12) %>%
  collapse_rows(columns = 1:4, valign = "top") 


##### Visualize the TEs on the chromosome level using a bar plot
#remove unclassified numbers 
REPEATS <- REPEATS[REPEATS$order!="Unclassified repeat",]

bar <- REPEATS %>%
  group_by(chr, order) %>%
  summarise(num = n()) %>%
  mutate(order == factor(order,
                         levels = c("LTR-RT", "LINE", "TIR", "MITE", "Helitron"
                         ))) %>%
  arrange(chr, order)

#plot 
colors <- c("green3","red","blue","cyan","magenta")

ggplot(bar,
       aes(x=fct_reorder(chr,parse_number(chr),.desc = F), 
           y=num, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors)+
  ggtitle("Number of transposable element identified in each chromosome - Hap2") + 
  xlab("Chromosomes") + ylab("TE count") + labs(fill="TE families")

#plot pie chart   
pie_genome_perc_hap1  <- 
  pie_genome_hap1 %>% mutate(perc=paste0(round(values/sum(values)*100,2),"%"," ","-"," ",pie_genome_hap1$order))

colors <- c("green3","red","blue","cyan","magenta","white")

pie(pie_genome_perc_hap1$values, labels = pie_genome_perc_hap1$perc, 
    main = "Proportion of transposable element identified occupied in the genome - TME117 Hap2",col = colors)
legend("topright", c("Helitron","LINE","LTR-RT","MIT","TIR","Unmasked"), cex = 0.8,fill = colors)


##### For Circos plot
hap2_repeats <- REPEATS
hap1_repeats <- REPEATS

################### Circos plot commands 
pacman::p_load(circlize,data.table)
setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/TME117-assembly/FinalAssembly/anno_Jan23/")

#compare hap1 and hap2 gene annotation
gff1 <- fread("Manihot_esculenta_Hap1.gff3", h = F)
gff2 <- fread("Manihot_esculenta_Hap2.gff3", h = F)

#extact chromosome, annotation names, start and end columns 
gff1 <- gff1[, .(V1, V4,V5,V3)]
colnames(gff1) <- c("chr","start","end", "name")
gff2 <- gff2[, .(V1, V4,V5,V3)]
colnames(gff2) <- c("chr","start","end", "name")

#length of the chromosomes
chrlen <- fread("../chromlen.hap2.txt")
colnames(chrlen) <- c("chr","len")
#get genes
genes1 <- gff1[name == "gene"]
genes2 <- gff2[name == "gene"]

# initialising the circle plot
# plot circos
circos.clear()
col_text <- "grey40"
circos.par("track.height" = 0.8, gap.degree = 1,cell.padding = c(0, 0, 0, 0))

circos.initialize(
  factors = c(
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18"
  ),
  xlim = matrix(c(rep(0, 18), chrlen$len), ncol = 2)
)

# chromosomes

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr,
              cex = 0.8, col = "black",
              facing = "outside", niceFacing = TRUE
  )
}, bg.col = "grey90", bg.border = F, track.height = 0.06)


circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h = "top", labels.cex = 0.5, col = col_text)
})

#plot
#####
circos.genomicDensity(hap1_repeats,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1000000,
                      col = c("#CC79A7"),
                      track.height = 0.08
) 

circos.genomicDensity(genes1,
                      type = "l",
                      lwd = 0.6,
                      count_by = "percent",
                      window.size = 1000000,
                      col = c("#D55E00"),
                      track.height = 0.08
)
