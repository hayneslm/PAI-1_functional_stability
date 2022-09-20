install.packages("tidyverse");
install.packages("reshape2");
library(tidyverse);
library(scales);
library(reshape2);

# Import  DESeq2 data
amplicons = c(1:12) 
data_48h = NULL
for (item in amplicons) {
  temp = read_tsv(paste('DeSeq48h/Amp_',item,'_Input_48Hr.txt',sep = '')) %>%
    mutate(Amplicon = item);
  data_48h = rbind(data_48h, temp)
}

# Reformat the data to make it more usable and to parse it out
data_48h = data_48h %>%separate(mutation, into = c("AA","Rest"),
                                sep = "(?<=[A-Z])(?=[0-9])") %>% 
  separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(mutation = paste(AA,Pos,Mut, sep = "")) %>%
  mutate(Pos = as.numeric(Pos)) %>% mutate(log2FoldChange = -log2FoldChange)

WT_48h = data_48h %>% filter(AA==Mut) %>% mutate(type = "wildtype")
data_48h_filtered=data_48h %>% na.omit() %>% anti_join(.,WT_48h,by = 'mutation')
nonsense_48h = data_48h_filtered %>% filter(Mut == 'X') %>% 
  mutate(type = "nonsense")
missense_48h = data_48h_filtered %>% filter(AA != Mut) %>%
  anti_join(.,nonsense_48h,by='mutation') %>% 
  mutate(type = "missense")
labeled_48h = rbind(WT_48h,missense_48h,nonsense_48h) %>% 
  mutate(padj = ifelse(padj>=0.05,"fail","pass" )) #pass/fail qvalue to add shapes to MA plot


# Generate an MA plot...setting basemean cutoff to 10.
ggplot() + 
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  ggtitle('48h MA Plot')+
  xlab("Base Mean Score")+
  ylab("Log2 Fold Change")+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(labeled_48h, type == "missense"), 
              aes(baseMean, log2FoldChange, color = type, shape = padj), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(labeled_48h, type == "nonsense"), 
              aes(baseMean, log2FoldChange, color = type, shape = padj), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(labeled_48h, type == "wildtype"),
            aes(baseMean, log2FoldChange, color = type, shape = padj),
           stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 10, color = "black", linetype = "dashed")+
  theme_classic()

# Filter based on cutoffs and genrate a HM of 48h enriched/depleted
filtered_missense = missense_48h %>% filter(padj < 0.05, baseMean > 10)

# adjust filtered missense to take into account WT selection
position = c(1:379)

temp = NULL
missense_WTadj = NULL
for (pos in position) {
  log2.WT = WT_48h %>% filter(Pos == pos) %>% select(log2FoldChange) %>% 
    as.numeric() 
  temp = filtered_missense %>% filter(Pos == pos) %>% 
    mutate(Log2wt = log2FoldChange-log2.WT)
  missense_WTadj = rbind(missense_WTadj, temp)
}
# rename the WT corrected log2FoldChange
missense_WTadj = missense_WTadj %>% select(-log2FoldChange) %>% 
  rename(log2FoldChange = Log2wt)

aa20 =  "AVLIMFYWRKHDESTNQGCP"
aa20 = unlist(strsplit(aa20, split = ""))

PAI1seq = "VHHPPSYVAHLASDFGVRVFQQVAQASKDRNVVFSPYGVASVLAMLQLTTGGETQQQIQAAMGFKIDDKGMAPALRHLYKELMGPWNKDEISTTDAIFVQRDLKLVQGFMPHFFRLFRSTVKQVDFSEVERARFIINDWVKTHTKGMISNLLGKGAVDQLTRLVLVNALYFNGQWKTPFPDSSTHRRLFHKSDGSTVSVPMMAQTNKFNYTEFTTPDGHYYDILELPYHGDTLSMFIAAPYEKEVPLSALTNILSAQLISHWKGNMTRLPRLLVLPKFSLETEVDLRKPLENLGMTDMFRQFQADFTSLSDQEPLHVAQALQKVKIEVNESGTVASSSTAVIVSARMAPEEIIMDRPFLFVVRHNPTGTVLFMGQVMEP"
Pos = c(1:nchar(PAI1seq))
AA = unlist(strsplit(PAI1seq, split=""))
WT_coded = cbind(Pos, AA, replicate(nchar(PAI1seq),10)) %>% as.data.frame() %>% 
  rename(value = V3) %>% mutate(value = as.numeric(value)) %>%
  rename(Mut = AA) %>% 
  mutate(padj = NA) %>% 
  as_tibble()

HM_data = missense_WTadj %>% 
  select(Pos, Mut, log2FoldChange,padj) %>% rename(value = log2FoldChange) 

# Filter out data that was depleted at t=0 and make Heat Map...I need to work on this
depleted_0h = read_tsv("all_data_2.txt") %>% 
  filter(log2FoldChange<=0, padj<0.05, baseMean > 10)

HM_functional = anti_join(HM_data,depleted_0h %>% 
                            rename(value = log2FoldChange), 
                          by = c("Pos","Mut"))

HM_dataF = HM_functional %>% select(Pos, Mut, value,padj) %>% 
  rbind(.,WT_coded)

all_data_sum = read_csv('alldata_sum.csv') %>% filter(Mut != 'X') %>% 
  select(mutation, summed) %>%
  separate(mutation, into = c("AA","Rest"),
           sep = "(?<=[A-Z])(?=[0-9])") %>% 
  separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(mutation = paste(AA,Pos,Mut, sep = "")) %>%
  mutate(Pos = as.numeric(Pos))

not_present = anti_join(all_data_sum, missense_WTadj, by = 'mutation') %>%
  mutate(value = NA) %>% 
  select(Pos, Mut, value) %>% 
  mutate(padj = NA)

HMlist = rbind(HM_dataF, not_present) %>% mutate(Pos = as.numeric(Pos)) %>% 
  filter(Pos<=379)

HMlist  <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)
levels(HMlist$Pos)


ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
                       values = rescale(c(-3, 0, 5, 10), 
                                        to = c(0,1)),
                       na.value = "ivory2") +
  theme(panel.background = element_rect(fill = 'white'),
        #legend.position = "none",
        axis.text=element_text(size=24))+
  ggtitle("Heat Map of Enriched/Depleted at 48h")

ggsave("HM_log2_legend.png", device = "png", width = 44, height = 8.5, units = "in")


#Next bit is to take the sum of values at each position
#First Melt the data and recast
HMfunctional.melt = melt(HM_functional,measure.vars = "value")
HMfunctional.cast = dcast(HMfunctional.melt,Mut~Pos) #%>% 
#mutate_all(~replace(., is.na(.), 0))
list.pos = HM_functional$Pos %>% unique() %>% as.tibble() %>% 
  rename(Pos = value)

all.positions = c(1:379) %>% as.tibble() %>% mutate(Sum = 0) %>% 
  rename(Pos = value)

sums48h_functional = colSums(Filter(is.numeric, HMfunctional.cast), na.rm = T) %>% as.tibble() %>% cbind(list.pos,.) %>% rename(Sum = value) %>%
  rbind(.,anti_join(all.positions,., by = "Pos"))

#export as a .txt file
write_tsv(sums48h_functional,"sums48h_functional.txt")

