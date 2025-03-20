setwd("~")
dat0 <- read.table("Lipids/relevant_ions_genes_v2.csv")
ids <- readRDS("Lipids/spot_slices_ids.rds")



#dat_c <- cor(t(dat0), method = "spearman")^3
#dat_c <- dat_c[upper.tri(dat_c)]

create_cor_matrix <- function(x, all=F, leaves=F,p_input =.1){
  
  if(isTRUE(all)){
    dat <- dat0
  } else{
    dat <- dat0[,ids==x]
  }
  
dat_c <- cor(t(dat), method = "spearman")^3

dat <- data.frame(row=rownames(dat_c)[row(dat_c)[upper.tri(dat_c)]], 
           col=colnames(dat_c)[col(dat_c)[upper.tri(dat_c)]], 
           corr=dat_c[upper.tri(dat_c)])
dat$idx1 <- str_extract(dat$row,"^p_|^n_"); dat$idx1[is.na(dat$idx1)] <- "g_"
dat$idx2 <- str_extract(dat$col,"^p_|^n_"); dat$idx2[is.na(dat$idx2)] <- "g_"

dat_combined <- dat 
colnames(dat_combined) <- c("from","to","width","idx1","idx2")
dat_combined$value <- dat_combined$width
dat_combined$color <- "purple"
dat_combined$color[dat_combined$width>0] <- "orange"

nodes <- data.frame(id = unique(sort(c(dat_combined$from,dat_combined$to))))
nodes$label <- nodes$id
nodes$color <- "green"
nodes$color[grepl("^n_", nodes$id)] <- "blue"
nodes$color[grepl("^p_", nodes$id)] <- "red"
nodes$border <- "black"
nodes$background <- "orange"
nodes$shadow <- T
nodes$font.color <- "white"
nodes$shape <- "box"

###

dat_combined <- dat_combined %>%
  filter(abs(width) > p_input) %>% filter(idx1 != idx2)  ### INITIAL FILTERING CONDITION
all_feat <- c(dat_combined$from, dat_combined$to)

nodes$size <- sapply(nodes$id, function(x) sum(all_feat==x))

nodes$font.size = 14
nodes$font.size[nodes$size>=10] = 16
nodes$font.size[nodes$size>=20] = 20

if(!isTRUE(leaves)){   ##### Focus on centroid hubs, or focus on highly connected network
  
nodes <- nodes %>% filter(size>2)
dat_combined <- dat_combined %>% filter(to %in% nodes$id , from %in% nodes$id,)
all_feat <- c(dat_combined$from, dat_combined$to)
nodes$size <- sapply(nodes$id, function(x) sum(all_feat==x))

} else{
centroids <- nodes %>% filter(size>2) %>% pull(id)
dat_combined <- dat_combined %>% filter((to %in% centroids) | (from %in% centroids))
all_feat <- c(dat_combined$from, dat_combined$to)
nodes$size <- sapply(nodes$id, function(x) sum(all_feat==x))
}

nodes <- nodes %>% filter(id %in% c(dat_combined$from,dat_combined$to)) %>% arrange(size)
nodes$size = 5*max(-log(exp(nodes$size)/max(nodes$size)))-(-log(exp(nodes$size)/max(nodes$size)))
nodes$value <- nodes$size
nodes$font.size <- round(nodes$size)
nodes$label <- str_remove(nodes$label,"\\.[0-9]+\\.[0-9]+$")
nodes$label <- str_remove(nodes$label,"^[a-z]{1}_")

net <- list(nodes=nodes, edges = dat_combined)

return(net)

}
unique(ids)

wt22 <- create_cor_matrix("WT_22mo",F)
wt12 <- create_cor_matrix("WT_12mo",F,p_input = .075)
wt6 <- create_cor_matrix("WT_6mo",F,F,p_input=.3)

tg22 <- create_cor_matrix("Tg_22mo",F,p_input=.275)
tg12 <- create_cor_matrix("Tg_12mo",F,p_input=.125)
tg6 <- create_cor_matrix("Tg_6mo",F,p_input=.10)
apoe3 <- create_cor_matrix("APOE3_12mo",F,p_input=.30)
apoe4 <- create_cor_matrix("APOE4_12mo",F,p_input=.15)

net0 <- create_cor_matrix("A",all = T,leaves=T,p_input = .1)


net <- wt12

visNetwork(net$nodes,net$edges, directed=T,width=2400,height=1600) %>%
  #  visLayout(hierarchical = T) %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
#  visNodes(scaling = list(label = list(enabled = T))) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
  #  visInteraction(navigationButtons = TRUE, 
  #                 dragNodes = T,
  #                 dragView = T, zoomView = T) %>%
  visInteraction(dragNodes = T, 
                 dragView = T, 
                 zoomView = T) %>%  
  visPhysics(solver = "repulsion", 
             repulsion = list("centralGravity"=.6,
                              "springLength"=300))%>%
  visLayout(randomSeed = 10) 


############

observed <- wt6$edges %>% mutate(pair = paste0(if_else(width>0,"P","N"),"_",to,"-",from)) %>% pull(pair)

if(isTRUE(all)){
  dat <- dat0
} else{
  dat <- dat0[,ids==x]
}

dat_c <- cor(t(dat), method = "spearman")^3
dat <- data.frame(row=rownames(dat_c)[row(dat_c)[upper.tri(dat_c)]], 
                  col=colnames(dat_c)[col(dat_c)[upper.tri(dat_c)]], 
                  corr=dat_c[upper.tri(dat_c)])
dat$idx1 <- str_extract(dat$row,"^p_|^n_"); dat$idx1[is.na(dat$idx1)] <- "g_"
dat$idx2 <- str_extract(dat$col,"^p_|^n_"); dat$idx2[is.na(dat$idx2)] <- "g_"
colnames(dat) <- c("from","to","width","idx1","idx2")

set.seed(99)
edges_fn <- function(x){

  if(x%%60==0){
    print(x)
  } 
dat$width <- sample(dat$width, nrow(dat),replace=F)

dat_combined <- dat 

nodes <- data.frame(id = unique(c(dat_combined$from,dat_combined$to)))

dat_combined <- dat_combined %>%
  filter(abs(width) > .1) %>% filter(idx1 != idx2)  ### INITIAL FILTERING CONDITION
all_feat <- c(dat_combined$from, dat_combined$to)

nodes <- nodes %>% left_join(data.frame(table(all_feat)) %>% rename("size"="Freq"),by=c("id"="all_feat"))

  nodes <- nodes %>% filter(size>2)
  dat_combined <- dat_combined %>% filter(to %in% nodes$id , from %in% nodes$id,)
  expected <- dat_combined %>% mutate(pair = paste0(if_else(width>0,"P","N"),"_",to,"-",from)) %>% pull(pair)
  return(sum(expected %in% observed))
#  centroids <- nodes %>% filter(size>2) %>% pull(id)
#  dat_combined <- dat_combined %>% filter((to %in% centroids) | (from %in% centroids))

}
nulls <- sapply(1:1000,edges_fn)


require(vcd)
require(MASS)

# data generation
ex <- rexp(10000, rate = 1.85) # generate some exponential distribution
control <- abs(rnorm(10000)) # generate some other distribution

# estimate the parameters
fit1 <- fitdistr(nulls, "normal") 
fit2 <- fitdistr(control, "exponential")

# goodness of fit test
ks.test(nulls, "pexp", fit1$estimate) # p-value > 0.05 -> distribution not refused
ks.test(control, "pexp", fit2$estimate) #  significant p-value -> distribution refused

# plot a graph
hist(nulls,freq=F)
curve(dnorm(x, mean = fit1$estimate[1],sd = fit1$estimate[2]), from = 0, col = "red", add = TRUE)
#############3

create_dat <- function(x){
dat <- dat0[,ids==x]

dat_c <- cor(t(dat), method = "spearman")^3
dat <- data.frame(row=rownames(dat_c)[row(dat_c)[upper.tri(dat_c)]], 
                  col=colnames(dat_c)[col(dat_c)[upper.tri(dat_c)]], 
                  corr=dat_c[upper.tri(dat_c)])
dat$idx1 <- str_extract(dat$row,"^p_|^n_"); dat$idx1[is.na(dat$idx1)] <- "g_"
dat$idx2 <- str_extract(dat$col,"^p_|^n_"); dat$idx2[is.na(dat$idx2)] <- "g_"
colnames(dat) <- c("from","to","width","idx1","idx2")
return(dat)
}

apoe3_dat <- create_dat("APOE3_12mo")
apoe4_dat <- create_dat("APOE4_12mo")

wt12_dat <- create_dat("WT_12mo")
tg12_dat <- create_dat("Tg_12mo")

perm <- function(dat,p_input){
  dat$width <- sample(dat$width, nrow(dat),replace=F)
  
  dat_combined <- dat 
  
  nodes <- data.frame(id = unique(c(dat_combined$from,dat_combined$to)))
  
  dat_combined <- dat_combined %>%
    filter(abs(width) > .1) %>% filter(idx1 != idx2)  ### INITIAL FILTERING CONDITION
  all_feat <- c(dat_combined$from, dat_combined$to)
  
  nodes <- nodes %>% left_join(data.frame(table(all_feat)) %>% rename("size"="Freq"),by=c("id"="all_feat"))
  
  nodes <- nodes %>% filter(size>2)
  dat_combined <- dat_combined %>% filter(to %in% nodes$id , from %in% nodes$id,)
  expected <- dat_combined %>% mutate(pair = paste0(if_else(width>0,"P","N"),"_",to,"-",from)) %>% pull(pair)
}


set.seed(99)
edges_fn_pair <- function(x,dat1,dat2,p1,p2){
  
  if(x%%60==0){
    print(x)
  } 

  perm1 <- perm(dat1, p1)
  perm2 <- perm(dat2,p2)
  
  stat_df <- data.frame(i = x, 
                        p1 = length(perm1),
                        p2 = length(perm2),
                        overlap = length(intersect(perm1,perm2))) %>%
    mutate(overlap.pct = overlap/(p1+p2-overlap))
  return(stat_df)

}
"FDE725FF"
nulls <- map_dfr(1:500,function(x) edges_fn_pair(x,wt12_dat,tg12_dat,.045,.125))
nulls <- map_dfr(1:500,function(x) edges_fn_pair(x,tg12_dat,apoe4_dat,.125,.15))

nulls %>% write.csv("APOE3_APOE4_null_dist.csv")
nulls %>% write.csv("WT12_TG12_null_dist.csv")

o1 <- tg12$edges %>% mutate(pair = paste0(if_else(width>0,"P","N"),"_",to,"-",from)) %>% pull(pair)
o2 <- apoe4$edges %>% mutate(pair = paste0(if_else(width>0,"P","N"),"_",to,"-",from)) %>% pull(pair)
o <- length(intersect(o1,o2))/(length(o1)+length(o2)+length(intersect(o1,o2)))
quantile(nulls$overlap.pct, by= seq(0,1,0.05))
o

plt <- nulls %>%
  ggplot(aes(x= overlap.pct))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(xintercept = o, linetype="dashed",color ="blue")+
  theme_bw()
plt
ggsave("TG12_APOE4_null_dist_network_edges_overlap.png", plt, height=6,width=5)

# wt22 <- create_cor_matrix("WT_22mo",F)
# wt12 <- create_cor_matrix("WT_12mo",F,p_input=.065)
# wt6 <- create_cor_matrix("WT_6mo",F,F,p_input=.15)
# 
# tg22 <- create_cor_matrix("Tg_22mo",F,p_input=.275)
# tg12 <- create_cor_matrix("Tg_12mo",F,p_input=.125)
# tg6 <- create_cor_matrix("Tg_6mo",F,p_input=.10)
# apoe3 <- create_cor_matrix("APOE3_12mo",F,p_input=.30)
# apoe4 <- create_cor_matrix("APOE4_12mo",F,p_input=.15)


library(tidyverse)
create_dat <- function(x){
  dat <- dat0[,ids==x]
  
  dat_c <- cor(t(dat), method = "spearman")^3
  dat <- data.frame(row=rownames(dat_c)[row(dat_c)[upper.tri(dat_c)]], 
                    col=colnames(dat_c)[col(dat_c)[upper.tri(dat_c)]], 
                    corr=dat_c[upper.tri(dat_c)])
  dat$idx1 <- str_extract(dat$row,"^p_|^n_"); dat$idx1[is.na(dat$idx1)] <- "g_"
  dat$idx2 <- str_extract(dat$col,"^p_|^n_"); dat$idx2[is.na(dat$idx2)] <- "g_"
  colnames(dat) <- c("from","to","width","idx1","idx2")
  return(dat)
}

apoe3_dat <- create_dat("APOE3_12mo")
apoe4_dat <- create_dat("APOE4_12mo")

wt12_dat <- create_dat("WT_12mo") %>% mutate(condition = "WT", age = "12mo", group = "WT_12mo") %>% filter(idx1!=idx2)
tg12_dat <- create_dat("Tg_12mo")%>% mutate(condition = "Tg", age = "12mo", group = "TG_12mo")%>% filter(idx1!=idx2)
wt6_dat <- create_dat("WT_6mo")%>% mutate(condition = "WT", age = "06mo", group = "WT_06mo")%>% filter(idx1!=idx2)
tg6_dat <- create_dat("Tg_6mo")%>% mutate(condition = "Tg", age = "06mo", group = "TG_06mo")%>% filter(idx1!=idx2)
wt22_dat <- create_dat("WT_22mo")%>% mutate(condition = "WT", age = "22mo", group = "WT_22mo")%>% filter(idx1!=idx2)
tg22_dat <- create_dat("Tg_22mo")%>% mutate(condition = "Tg", age = "22mo", group = "TG_22mo")%>% filter(idx1!=idx2)

rbind(wt12_dat,tg12_dat,wt6_dat,tg6_dat,wt22_dat,tg22_dat) %>%
  ggplot(aes(x = width, color =group))+
  geom_histogram(fill="white")+
  theme_bw()+
  facet_grid(age~condition,scales="free")

all_dat <- cbind(wt12_dat[,1:2], wt12_dat$width, tg12_dat$width,wt6_dat$width,tg6_dat$width,
                 wt22_dat$width,tg22_dat$width)
all_dat <- all_dat %>% pivot_longer(-any_of(c("from","to")), names_to="group", values_to="correlation") %>%
  mutate(group = str_remove(group,"\\$width")) %>%
  mutate(condition = str_to_upper(str_extract(group,"[a-z]+"))) %>%
  mutate(age = paste0(str_extract(group,"[0-9]+"),"mo")) %>%
  mutate(age = factor(age, levels = c("6mo","12mo","22mo")))

median_lists <- all_dat %>%
  group_by(group) %>%
  summarize(middle = median(correlation,na.rm=T)) %>%
  ungroup()


all_dat1 <- all_dat %>% mutate(id = paste0(from,"-",to)) %>% 
  filter(!is.na(correlation)) %>%
  group_by(id,condition) %>%
  mutate(count = dplyr::n()) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(values = sum(!is.na(correlation))) %>%
#  mutate(rank = order(correlation)/values) %>%
  mutate(rank = correlation) %>%
  ungroup()
wt_all3 <- all_dat1 %>% filter(count==3,condition=="WT") %>% pull(id) %>% unique()
tg_all3 <- all_dat1 %>% filter(count==3,condition=="TG") %>% pull(id) %>% unique()

#slopes <- all_dat1 %>% group_by(condition, from,to) %>% summarize(rank = (max(rank)-min(rank))/3)
slopes <- all_dat1 %>% 
  group_by(condition,id) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(r = cor(rank,age)) %>%
 # do(
    #lm = lm(rank~age,data=.)
##    r = cor(rank,age)
#  ) %>%
  ungroup()

saveRDS(slopes, file ="Lipids/Ranked_Correlations_Age_Ion_Gene_Interactions.rds")
slopes <- readRDS(file ="Lipids/Raw_Correlations_Age_Ion_Gene_Interactions.rds")

#slopes %>% write.csv("Lipids/Raw_Correlations_Age_Ion_Gene_Interactions.csv")

ids <- slopes %>% mutate(id = paste0(from,"-",to)) %>% filter(abs(r) < .999, id %in% intersect(wt_all3,tg_all3)) %>%
  top_n(108, abs(r)) %>% pull(id) %>% unique()
length(ids)


ids <- slopes %>% mutate(id = paste0(from,"-",to)) %>% filter(abs(r) < .999, id %in% intersect(wt_all3,tg_all3)) %>%
  filter(abs(correlation) >= .01) %>% 
  top_n(10000, abs(r)) %>% arrange(desc(abs(r))) %>% pull(id) %>% unique()

res <- all_dat1 %>%
  filter(id %in% ids) %>%
  group_by(id) %>%
  mutate(age = as.numeric(age)) %>%
  do(
    nullm = lm(rank~age*condition,data=.),
    fullm = lm(rank ~ age,data=.),
  ) %>%
  ungroup()
res$p <- sapply(1:nrow(res), function(x) anova(res$nullm[[x]],res$fullm[[x]])$`Pr(>F)`[2])
res$p <- p.adjust(res$p, "BH")
summary(res$p)
p_ids <- ids[ids %in% res$id[res$p <= .1]] %>% head(36)

plt <- all_dat1 %>%
  filter(id %in% p_ids) %>%
  mutate(count = factor(count, levels=c(2,3), labels = c("2","3"))) %>%
  mutate(alph = ifelse(id %in% ids, .8,.2)) %>%
  mutate(ids2 = ifelse(id %in% ids, id, NA)) %>%
  ggplot(aes(x=age, y= rank,alpha=alph, group = condition,shape=condition,color =condition, linetype=count))+
  geom_point(size=4, alpha=1)+
  geom_line(alpha=1)+
  theme_bw()+
  ylab("Correlation")+
  scale_color_igv()+
  scale_linetype_manual(values = c("2"="dashed","3"="solid"))+
  theme(axis.text = element_text(color = "black"))+
  guides(alpha = "none", color = "none")+
  facet_wrap(from~to, ncol=6,scales="free_y")#+
  #labs(caption = "Y-axis are normalized ranks of sorted ion-gene correlations for each pallidum.")
plt
ggsave("TOP_36_correlated_interactions_raw_corr_values_v3.png",plt, height=14,width=15)

library(ggsci)



#############
setwd("~")
library(ggsci)
assignments <- read.csv("module_assignments_121323.csv")
eigen_mat <- read.csv("module_eigenvalues_121323.csv")
ids <- readRDS("Lipids/spot_slices_ids.rds")

cor2 <- function(...){
  
  cor(...,method="spearman")
}

eigen_mat$id <- ids
eigen_cor <- eigen_mat %>%
  select(-X) %>%
  relocate(id) %>%
  split(.$id) %>% 
  map(select, -c(id)) %>% 
  map(cor2)

eigen_cor <- map_dfr(1:8, function(x){
  dat_c <- eigen_cor[[x]]
  dat <- data.frame(row=rownames(dat_c)[row(dat_c)[upper.tri(dat_c)]], 
                    col=colnames(dat_c)[col(dat_c)[upper.tri(dat_c)]], 
                    corr=dat_c[upper.tri(dat_c)])
  dat$id <- names(eigen_cor)[x]
  dat <- dat %>% relocate(id)
  return(dat)
}
  )
eigen_cor$V1 <- str_extract(eigen_cor$row, "^G\\.|^N\\.|^P\\.")
eigen_cor$V2 <- str_extract(eigen_cor$col, "^G\\.|^N\\.|^P\\.")

eigen_cor <- eigen_cor %>% filter(V1!=V2) %>%
  mutate(age = str_extract(id, "6mo|12mo|22mo")) %>%
  mutate(age = factor(age, levels = c("6mo","12mo","22mo"))) %>%
  mutate(condition = str_extract(id,"Tg|WT|APOE3|APOE4")) %>%
  group_by(row,col) %>%
  mutate(max_cor = max(abs(corr))) %>%
  ungroup()

plt <- eigen_cor %>%
  filter(!grepl("APOE",condition)) %>%
  filter(max_cor > .5) %>%
  filter(!grepl("grey", row), !grepl("grey",col)) %>%
ggplot(aes(x=age, y= corr, group = condition,shape=condition,color =condition))+
  geom_point(size=4, alpha=1)+
  geom_line(alpha=1)+
  theme_bw()+
  ylab("p")+
  scale_color_igv()+
  theme(axis.text = element_text(color = "black"))+
  guides(alpha = "none", color = "none")+
  facet_wrap(row~col, ncol=7,scales="free_y")#+
#labs(caption = "Y-axis are normalized ranks of sorted ion-gene correlations for each pallidum.")
#plt  
ggsave("Mixed_modalities_rho_.5.png", plt, height=14,width=19,dpi=1500)

library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")

color_df <- map_dfr(unique(assignments$module[!grepl("^p_|^n_",assignments$feature)]), function(y){
  print(y)
  GO_results8 <- enrichGO(gene = assignments %>% filter(!grepl("^p_|^n_", feature), module ==y) %>% pull(feature) ,
                          OrgDb = "org.Mm.eg.db", keyType = "SYMBOL",
                          ont = "all")
  if(!is.null(GO_results8)){
  sig_cat <- GO_results8@result %>% separate(BgRatio,sep="/",c("n_genes","n_back")) %>%
    mutate(n_genes = as.numeric(n_genes), 
           n_back = as.numeric(n_back)) %>%
    mutate(BgRatio = n_genes/n_back) %>% dplyr::select(-any_of(c("n_genes","n_back"))) %>%
 #   mutate(Description = str_wrap(Description,20)) %>%
    filter(Count >= 4) %>% mutate(Module = y)
  return(sig_cat)
  }
})


color_df %>%
  group_by(Module) %>% 
  ## keep first column only and name it 'keywords':
  select(Module,'keywords' = 3) %>%
  ## multiple cell values (as separated by a blank)
  ## into separate rows:
  separate_rows(keywords, sep = " ") %>%
  group_by(Module,keywords) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

GO_blue <- enrichGO(gene =nick_hits, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "all")

library(clusterProfiler)
library(org.Hs.eg.db)


#genes <- read.csv("sample_genes.csv")
GO_blue <- enrichGO(gene = assignments %>%
                          filter(!grepl("^p_|^n_", feature), module =="blue") %>%
                          pull(feature) %>% str_to_upper(), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "all")
GO_turq <- enrichGO(gene = assignments %>%
                      filter(!grepl("^p_|^n_", feature), module =="turquoise") %>%
                      pull(feature) %>% str_to_upper(), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "all")
GO_brown <- enrichGO(gene = assignments %>%
                      filter(!grepl("^p_|^n_", feature), module =="brown") %>%
                      pull(feature) %>% str_to_upper(), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "all")

blue_go <- GO_blue@result %>% separate(BgRatio,sep="/",c("n_genes","n_back")) %>%
  mutate(n_genes = as.numeric(n_genes), 
         n_back = as.numeric(n_back)) %>%
  mutate(BgRatio = n_genes/n_back) %>% dplyr::select(-any_of(c("n_genes","n_back"))) %>%
  mutate(Description = str_wrap(Description,20)) %>%
  filter(Count >= 4) 

turq_go <- GO_turq@result %>% separate(BgRatio,sep="/",c("n_genes","n_back")) %>%
  mutate(n_genes = as.numeric(n_genes), 
         n_back = as.numeric(n_back)) %>%
  mutate(BgRatio = n_genes/n_back) %>% dplyr::select(-any_of(c("n_genes","n_back"))) %>%
  mutate(Description = str_wrap(Description,20)) %>%
  filter(Count >= 4) 

brown_go <- GO_brown@result %>% separate(BgRatio,sep="/",c("n_genes","n_back")) %>%
  mutate(n_genes = as.numeric(n_genes), 
         n_back = as.numeric(n_back)) %>%
  mutate(BgRatio = n_genes/n_back) %>% dplyr::select(-any_of(c("n_genes","n_back"))) %>%
  mutate(Description = str_wrap(Description,20)) %>%
  filter(Count >= 4) 


blue_go %>% write.csv("genes_blue_module_go_results.csv")
turq_go %>% write.csv("genes_turq_module_go_results.csv")
brown_go %>% write.csv("genes_brown_module_go_results.csv")

getwd()
head(eigen_cor)

create_cor_matrix <- function(x, all=F, leaves=F,p_input =.1){
  


  dat_combined <- eigen_cor %>% filter(id=="WT_12mo")
  dat_combined$corr <- dat_combined$corr^2
  colnames(dat_combined)[1:6] <- c("id","from","to","width","idx1","idx2")
  dat_combined$value <- dat_combined$width
  dat_combined$color <- "purple"
  dat_combined$color[dat_combined$width>0] <- "orange"
  
  nodes <- data.frame(id = unique(sort(c(dat_combined$from,dat_combined$to))))
  nodes$label <- nodes$id
  nodes$color <- "green"
  nodes$color[grepl("N", nodes$id)] <- "blue"
  nodes$color[grepl("P", nodes$id)] <- "red"
  nodes$border <- "black"
  nodes$background <- "orange"
  nodes$shadow <- T
  nodes$font.color <- "white"
  nodes$shape <- "box"
  
  ###
  
  dat_combined <- dat_combined %>%
    select(-id) %>%
    filter(abs(width) > .1) %>% filter(idx1 != idx2)  ### INITIAL FILTERING CONDITION
  all_feat <- c(dat_combined$from, dat_combined$to)
  
  nodes$size <- sapply(nodes$id, function(x) sum(all_feat==x))
  
  nodes$font.size = 14
  nodes$font.size[nodes$size>=10] = 16
  nodes$font.size[nodes$size>=20] = 20
  
 # if(!isTRUE(leaves)){   ##### Focus on centroid hubs, or focus on highly connected network
    
    nodes <- nodes %>% filter(size>2)
    dat_combined <- dat_combined %>% filter(to %in% nodes$id , from %in% nodes$id,)
    all_feat <- c(dat_combined$from, dat_combined$to)
    nodes$size <- sapply(nodes$id, function(x) sum(all_feat==x))
    
  # } else{
  #   centroids <- nodes %>% filter(size>2) %>% pull(id)
  #   dat_combined <- dat_combined %>% filter((to %in% centroids) | (from %in% centroids))
  #   all_feat <- c(dat_combined$from, dat_combined$to)
  #   nodes$size <- sapply(nodes$id, function(x) sum(all_feat==x))
  # }
  
  nodes <- nodes %>% filter(id %in% c(dat_combined$from,dat_combined$to)) %>% arrange(size)
  nodes$size = 5*max(-log(exp(nodes$size)/max(nodes$size)))-(-log(exp(nodes$size)/max(nodes$size)))
  nodes$value <- nodes$size
  nodes$font.size <- round(nodes$size)+2
 # nodes$label <- str_remove(nodes$label,"\\.[0-9]+\\.[0-9]+$")
#  nodes$label <- str_remove(nodes$label,"^[a-z]{1}_")
  
  net <- list(nodes=nodes, edges = dat_combined)
  
  visNetwork(net$nodes,net$edges, directed=T,width=800,height=700) %>%
    #  visLayout(hierarchical = T) %>%
    visIgraphLayout(layout = "layout_with_fr") %>%
    #  visNodes(scaling = list(label = list(enabled = T))) %>%
    visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
    #  visInteraction(navigationButtons = TRUE, 
    #                 dragNodes = T,
    #                 dragView = T, zoomView = T) %>%
    visInteraction(dragNodes = T, 
                   dragView = T, 
                   zoomView = T) %>%  
    visPhysics(solver = "repulsion", 
               repulsion = list("centralGravity"=.6,
                                "springLength"=300))%>%
    visLayout(randomSeed = 10) 
  return(net)
  
}
