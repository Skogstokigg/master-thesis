###Libraries
library("fpc")
require(maps)
library("vegan")
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggtext)
library(stats)


### 1. Metadata
metadata<-read.table(file="/mnt/lustre/scratch/elopez/0_Previous_Analysis/Selected_Malaspina_Profiles_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(metadata)<-metadata$MPCode
metadata$Sample<-metadata$MPCode


### 2. Sampling Site Map
col_pal<-brewer.pal(11,"Spectral")
col_grad<-rev(colorRampPalette(col_pal[4:9])(n=299))

#Load world map data
world_map <- map_data("world")


fig_map<-ggplot(world_map, aes(x = long, y = lat, group = group))+geom_polygon(fill="lightgray", colour = "lightgray")+geom_point(metadata,mapping=aes(y=Lat,x=Long,group=MPCode,col=Depth_m,size=Depth_m))+coord_cartesian(xlim=c(-180,180),ylim=c(-90,90))+geom_text(data=metadata,size=4,aes(y=Lat+5,x=Long,label=Station,group=Station))+theme_bw()+xlab("Longitude")+ylab("Latitude")+scale_colour_gradientn(colours =col_grad,limits = c(4000, 0),name="Depth",trans = 'reverse')+guides(size = "none")

ggsave("Malaspina_Profiles_Figure.pdf",plot=fig_map,device=cairo_pdf,width=14,height=7,pointsize=8)



### 3. Scaffolds Length Distribution
file1 <- read.table(file="/mnt/lustre/scratch/elopez/0_Previous_Analysis/CountsxLength/length_count.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
file1 <- aggregate(. ~ Length, data = file1, sum)
rownames(file1) <- file1$Length
file1$Length <- file1$Length

fig_length<- ggplot(file1, aes(x = Length)) +
  geom_histogram(bins = 70, color = "black", size = 0.6, fill = "coral") +
  theme(plot.le.x = element_text(angle = 90,hjust = 1,size=4),
        title = element_markdown(),
        subtitle = element_markdown()) +
  labs(title = "", y = "Count", x = "Scaffold Length") + 
  scale_x_continuous(labels = comma) +
  theme_bw()

ggsave("Scaffolds_Length_Distribution.pdf", plot = fig_length, device = cairo_pdf, width = 6, height = 4, pointsize = 15)



### 4. RDA
rda_metadata<-subset(metadata,select=c(Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,layer)) # selected variables
rda_metadata<-na.omit(rda_metadata)


sample_layer<-as.vector(rda_metadata$layer)
sample_layer[which(sample_layer == "Epi")]<-"Epipelagic"
sample_layer[which(sample_layer == "Meso")]<-"Mesopelagic"
sample_layer[which(sample_layer == "Bathy")]<-"Bathypelagic"

rda_metadata$layer<-NULL

sdata<-as.data.frame(scale(rda_metadata,center=TRUE,scale=TRUE))

rdadata<-rda(sdata)

xlabel<-round((summary(rdadata)$cont$importance[2,1]*100),digits=1)
xlabel<-paste("Axis 1 (",xlabel,"% explained)",sep="")

ylabel<-round((summary(rdadata)$cont$importance[2,2]*100),digits=1)
ylabel<-paste("Axis 2 (",ylabel,"% explained)",sep="")

pdf("Malaspina_Profiles_RDA.pdf",width=6,height=6,pointsize=8)
biplot(rdadata,cex=rep((2), 2),xlab=xlabel,ylab=ylabel)
ordihull(rdadata,group =sample_layer,col = c("#081D58","#FEB24C","#41B6C4"))
legend('topright', legend=c("Epipelagic","Mesopelagic","Bathypelagic"), col=c("#FEB24C","#41B6C4","#081D58"), pch = 16)		 
dev.off()



### 5. Virus NMDS (beta-diversity)
vir_scaff_abund<-read.table(file="/mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

fdata<-vir_scaff_abund
rownames(fdata)<-fdata[,1]
colnames(fdata)[which(colnames(fdata) == "MP2233")]<-"MP2233bis"
fdata<-fdata[,-1]

dist_metric<-"bray"

# ---------------------------- all samples -------------------------------------------------------------------
dists<-vegdist(t(fdata), method = dist_metric) # calculation of Bray-Curtis dissimilarity
set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999) # perform the NMDS analysis to reduce dimensions
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-merge(data.scores,metadata,by="Sample",all.x=TRUE)

# Definition and colour of the water layers
mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))
zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")


# represent the result in a dotplot
fig_nmds<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(shape=18,size=2,alpha=0.90,aes(color=Zone))+theme_bw()+theme(text=element_text(size=16))+scale_colour_manual(name="Zone",values=zone_coloring)

ggsave("Malaspina_Profiles_Viruses_NMDS.pdf",plot=fig_nmds,device=cairo_pdf,width=9,height=7,pointsize=8)

# Correlogram - Pearson
library("GGally")

cgram_metadata<-subset(mdata,select=c(NMDS1,NMDS2,Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,All_BT,all_virus))

pdf("Malaspina_Profiles_Metadata_Correlogram.pdf",width=25,height=25,pointsize=8)
plot<-ggpairs(cgram_metadata)
print(plot)
dev.off()


# Correlogram - Spearman
p1 = cgram_metadata %>%
  ggpairs(.,lower = list(continuous = wrap("points", size=0.6)),
          upper = list(continuous = wrap("cor", method = "spearman",size= 2.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), 
		axis.text.y = element_text(angle = 0, hjust = 1, size=6),
		strip.background = element_rect(color = "grey80", fill = "grey80"),
		strip.text = element_text(colour = "black", size = 6))
		  
ggsave("Malaspina_Profiles_Metadata_Correlogram_Spearman.pdf",plot=p1,device=cairo_pdf,width=9,height=9,pointsize=10)


# An NMDS for each layer. Also a correlogram.
# -------------------------------------------- Epipelagic ------------------------------------------------------

epi_metadata<-metadata[metadata$layer == "Epi", ] # filter only epipelagic samples from metadata 
library(dplyr)
epi_fdata<-select(fdata,epi_metadata$Sample) # filter only epipelagic samples from fdata (vir_scaff_abund) 

dists<-vegdist(t(epi_fdata), method = dist_metric)
set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-merge(data.scores,epi_metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic"))


col_pal<-brewer.pal(11,"Spectral")
col_grad<-rev(colorRampPalette(col_pal)(n=299))

fig_nmds<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(shape=18,size=4,alpha=0.75,aes(color=Depth_m))+theme_bw()+theme(text=element_text(size=16))+scale_colour_gradientn(colours =col_grad,trans="reverse")+labs(color="Depth (m)") # legend with the most significant variable from the correlogram

ggsave("Malaspina_Profiles_Viruses_NMDS_Epi.pdf",plot=fig_nmds,device=cairo_pdf,width=9,height=7,pointsize=8)


# Correlogram
cgram_metadata_epi<-subset(mdata,select=c(NMDS1,NMDS2,Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,All_BT,all_virus))

pdf("Malaspina_Profiles_Metadata_Correlogram_Epi.pdf",width=25,height=25,pointsize=8)
plot<-ggpairs(cgram_metadata_epi)
print(plot)
dev.off()

# ------------------------------------------- Mesopelagic ------------------------------------------------------
meso_metadata<-metadata[metadata$layer == "Meso", ] # filter only epipelagic samples from metadata 
meso_fdata<-select(fdata,meso_metadata$Sample) # filter only epipelagic samples from fdata (vir_scaff_abund) 

dists<-vegdist(t(meso_fdata), method = dist_metric)
set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-merge(data.scores,meso_metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Mesopelagic"))

col_pal<-brewer.pal(11,"Spectral")
col_grad<-rev(colorRampPalette(col_pal)(n=299))	

fig_nmds<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+
			geom_point(shape=18,size=4,alpha=0.75,aes(color=NO3_WOA13))+
			theme_bw()+theme(text=element_text(size=16))+
			scale_colour_gradientn(colours =col_grad,trans="reverse")+
			labs(color="NO3 (µmol/L)") 

ggsave("Malaspina_Profiles_Viruses_NMDS_Meso.pdf",plot=fig_nmds,device=cairo_pdf,width=9,height=7,pointsize=8)


# Correlogram
cgram_metadata_meso<-subset(mdata,select=c(NMDS1,NMDS2,Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,All_BT,all_virus))

pdf("Malaspina_Profiles_Metadata_Correlogram_Meso.pdf",width=25,height=25,pointsize=8)
plot<-ggpairs(cgram_metadata_meso)
print(plot)
dev.off()



# ------------------------------------------- Bathypelagic -----------------------------------------------------
bathy_metadata<-metadata[metadata$layer == "Bathy", ] # filter only epipelagic samples from metadata 
bathy_fdata<-select(fdata,bathy_metadata$Sample) # filter only epipelagic samples from fdata (vir_scaff_abund) 

dists<-vegdist(t(bathy_fdata), method = dist_metric)
set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-merge(data.scores,bathy_metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Bathypelagic"))

col_pal<-brewer.pal(11,"Spectral")
col_grad<-rev(colorRampPalette(col_pal)(n=299))

fig_nmds<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+
			geom_point(shape=18,size=4,alpha=0.75,aes(color=SiO4_WOA13))+
			theme_bw()+theme(text=element_text(size=16))+
			scale_colour_gradientn(colours =col_grad,trans="reverse")+
			labs(color="SiO4 (µmol/L)") 

ggsave("Malaspina_Profiles_Viruses_NMDS_Bathy.pdf",plot=fig_nmds,device=cairo_pdf,width=9,height=7,pointsize=8)


# Correlogram
cgram_metadata_bathy<-subset(mdata,select=c(NMDS1,NMDS2,Depth_m,Temp,Sal_PSU,NO3_WOA13,PO4_WOA13,SiO4_WOA13,Oxygen,All_BT,all_virus))

pdf("Malaspina_Profiles_Metadata_Correlogram_Bathy.pdf",width=25,height=25,pointsize=8)
plot<-ggpairs(cgram_metadata_bathy)
print(plot)
dev.off()



####ANOSIM
vir_scaff_abund<-read.table(file="/mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(vir_scaff_abund)<-vir_scaff_abund[,1]
vir_scaff_abund<-vir_scaff_abund[,-1]


set.seed(666)
result_anosim<-anosim(t(as.matrix(vir_scaff_abund)), metadata[colnames(vir_scaff_abund),]$layer, distance = "bray", permutations = 1000)
result_anosim





### 6. Shannon's index (alpha-diversity)
# Load file with RPKMs
vir_scaff_abund<-read.table(file="/mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/RPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

fdata<-vir_scaff_abund
rownames(fdata)<-fdata[,1]
colnames(fdata)[which(colnames(fdata) == "MP2233")]<-"MP2233bis"
fdata<-fdata[,-1]

# Calculate diversity 
ind_metric <- "shannon"
Shannons_H<- diversity(t(fdata), index = ind_metric) 
data.scores<-as.data.frame(scores(Shannons_H))
data.scores$Sample<-rownames(data.scores)

alpha_diversity<-as.data.frame(Shannons_H)
alpha_diversity_number_samples <- nrow(alpha_diversity)
alpha_diversity_names <- rownames(alpha_diversity)
alpha_diversity_data <- as.data.frame(cbind(alpha_diversity_names, alpha_diversity))

# Depth
mdata<-merge(data.scores,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"

mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

zone_coloring<-c(brewer.pal(9,"YlOrRd")[c(4)],brewer.pal(9,"YlGnBu")[c(5,9)])
names(zone_coloring)<-c("Epipelagic","Mesopelagic","Bathypelagic")

# Plot
fig_alpha <- ggplot(data=mdata,aes(x=Dim1,y=Depth_m))+
				  geom_point(shape=16,size=3,alpha=0.75,aes(color=Depth_m),show.legend = FALSE) +
				  labs(x="Alpha-diversity (Shannon's H)", y="Sample Depth (m)") +
				  scale_color_gradientn(colors = (colorRampPalette(c("#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"))(n=299))) +
				  scale_y_reverse() +
				  geom_smooth(orientation = "y") +
				  theme_classic()
 
ggsave("Malaspina_Profiles_Viruses_Shannon.pdf",plot=fig_alpha,device=cairo_pdf,width=4,height=7,pointsize=8)
  
  


### 7. Virus abundance
# -------------------------------------------- Host Phylum -----------------------------------------------------
vir_host_phylum_abund<-read.table(file="/mnt/lustre/scratch/elopez/7_normalization/host_filtered_FPKM/Filtered_by_Predicted_Host_Score-Predicted_Host_2_Phylum-FPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mdata<-melt(vir_host_phylum_abund)
colnames(mdata)<-c("Taxon","Sample","Abundance")
summary(mdata)

mdata<-merge(mdata,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"

mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata<-mdata

host_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8])
names(host_coloring)<-c("Actinobacteria","Alphaproteobacteria","Cyanobacteria","Euryarchaeota","Bacteroidetes","Deinococcus−Thermus","Verrucomicrobia","Firmicutes","Fusobacteria","Planctomycetes","Gammaproteobacteria","Betaproteobacteria","Crenarchaeota","Unknown")

host_coloring<-host_coloring[c(as.vector(unique(fmdata$Taxon)))]

fig_host<-ggplot(fmdata,aes(y=Abundance,x=Sample,fill=Taxon))+
		geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+
		coord_flip()+theme_bw()+ labs(y = "Relative Abundance (RPKM)")+
		theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+
		scale_fill_manual(name="Host Phylum",values=host_coloring)+facet_grid(Zone ~ ., scales="free_y")

ggsave("Malaspina_Profiles_Virus_RaFAH_Filtered_0.14_Host_Phylum_AbundancexSample.pdf",plot=fig_host,device=cairo_pdf,width=10,height=12,pointsize=8)

# -------------------------------------------- Virus Family -----------------------------------------------------
vir_fam_abund<-read.table(file="/mnt/lustre/scratch/elopez/7_normalization/virus_unfiltered_FPKM/Filtered_by_Predicted_Host_Score-family-FPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mdata<-melt(vir_fam_abund)
colnames(mdata)<-c("Taxon","Sample","Abundance")
summary(mdata)

mdata<-merge(mdata,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata<-mdata[which(mdata$Abundance >= 1000),]


vir_fam_coloring<-c(brewer.pal(11,"Spectral"),brewer.pal(8,"Dark2")[8])
names(vir_fam_coloring)<-c("Lavidaviridae","Myoviridae","Microviridae","Inoviridae","Baculoviridae","Globuloviridae","Phycodnaviridae","Sphaerolipoviridae","Podoviridae","Mimiviridae","Siphoviridae","Unknown")

vir_fam_coloring<-vir_fam_coloring[sort(as.vector(unique(fmdata$Taxon)))]


fig_family<-ggplot(fmdata,aes(y=Abundance,x=Sample,fill=Taxon))+
				geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+
				coord_flip()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+
				facet_grid(Zone ~ ., scales="free_y")+ labs(y = "Relative Abundance (RPKM)")+
				scale_fill_manual(name="Virus Family",values=vir_fam_coloring)

ggsave("Malaspina_Profiles_Virus_VPFClass_Family_AbundancexSample.pdf",plot=fig_family,device=cairo_pdf,width=10,height=12,pointsize=8)


#-------------------------------------------- Virus Genus -----------------------------------------------------
vir_genus_abund<-read.table(file="/mnt/lustre/scratch/elopez/7_normalization/virus_unfiltered_FPKM/Filtered_by_Predicted_Host_Score-genus-FPKM.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mdata<-melt(vir_genus_abund)
colnames(mdata)<-c("Taxon","Sample","Abundance")
summary(mdata)

mdata<-merge(mdata,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fmdata<-mdata[which(mdata$Abundance >= 100),]


col_pal<-brewer.pal(11,"Spectral")
col_grad<-rev(colorRampPalette(col_pal)(n=299))

fig_genus<-ggplot(fmdata,aes(fill=log10(Abundance),y=Sample,x=Taxon))+
				geom_raster()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+
				labs(x = "Virus Genus")+
				facet_grid(Zone ~ ., scales="free")+scale_fill_gradientn(colours =col_grad)

ggsave("Malaspina_Profiles_Virus_VPFClass_Genus_AbundancexSample.pdf",plot=fig_genus,device=cairo_pdf,width=25,height=15,pointsize=8)





### 8. mTAGs Phylum abundance
# INPUT FILES
# Phylum
all_phylum_abund<-read.table(file="/mnt/lustre/bio/shared/malaspina/cnag-metagenomes/malaspina.mTags/output/2019-09-30.1_malaspina-160-metagenomes-cnag-SILVA132-0.99-merged-tables/Merged_table_phylum.txt",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(all_phylum_abund)<-gsub("\\.(.)+$","",colnames(all_phylum_abund),perl=TRUE) # modify column names
rownames(all_phylum_abund)<-all_phylum_abund$X # modify row names
all_phylum_abund$X<-NULL

prok_phylum_abund<-all_phylum_abund[grepl("Bacteria|Archaea", rownames(all_phylum_abund)),] # create a new variable with all the information of Bacteria and Archaea
prok_phylum_abund<-as.data.frame(t((t(prok_phylum_abund)/colSums(prok_phylum_abund))*100)) # normalise abundance by the sample size

taxon<-data.frame(do.call('rbind', strsplit(as.character(rownames(prok_phylum_abund)),';',fixed=TRUE))) # remove ";" separator
rownames(prok_phylum_abund)<-taxon$X2 # select phylum as rowname
prok_phylum_abund<-t(prok_phylum_abund) # transpose matrix


# Merge files
mdata<-merge(prok_phylum_abund,subset(metadata,select=c("layer")),by="row.names",all.x=TRUE)
mdata<-melt(mdata,id=c("Row.names","layer"))
colnames(mdata)<-c("Sample","Zone","Taxon","Abundance")
mdata$Sample<-as.character(mdata$Sample)
fmdata<-mdata[which(mdata$Abundance >= 1 & mdata$Zone != "NA"),] # select samples from vertical profiles

# Set names and colours
host_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(9,"Set1")[2],brewer.pal(8,"Accent"))
names(host_coloring)<-c("Actinobacteria","Alphaproteobacteria","Cyanobacteria","Euryarchaeota","Bacteroidetes","Deinococcus−Thermus","Verrucomicrobia","Firmicutes","Fusobacteria","Planctomycetes","Gammaproteobacteria","Betaproteobacteria","Crenarchaeota","Unknown","Proteobacteria","Thaumarchaeota","Chloroflexi","Gemmatimonadetes","Nitrospinae","WPS-2","Patescibacteria","Acidobacteria","Chlamydiae")

for (taxon in unique(fmdata$Taxon)) {
	if (taxon %in% names(host_coloring)) {
	} else {
		print(paste(taxon," missing from color palette!",sep=""))
	}
}
host_coloring<-host_coloring[c(as.vector(unique(fmdata$Taxon)))]


fmdata$Zone<-gsub("Epi","Epipelagic",fmdata$Zone)
fmdata$Zone<-gsub("Meso","Mesopelagic",fmdata$Zone)
fmdata$Zone<-gsub("Bathy","Bathypelagic",fmdata$Zone)

fmdata$Zone<-factor(fmdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

fig_bact<-ggplot(fmdata,aes(y=Abundance,x=Sample,fill=Taxon))+
			geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+
			coord_flip()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size=9))+
			scale_fill_manual(name="Phylum",values=host_coloring)+facet_grid(Zone ~ ., scales="free_y")

ggsave("Malaspina_Profiles_mTAGs_Phylum_AbundancexSample.pdf",plot=fig_bact,device=cairo_pdf,width=10,height=12,pointsize=8)




### 9. KEGG Module Abundance Barplot
# -------------------------------------------- AMG metabolism ----------------------------------------------------
metab_abund<-read.table(file="/mnt/lustre/scratch/elopez/7_normalization/AMG_abundance/Metabolism_Abundance.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

fdata<-metab_abund
rownames(fdata)<-fdata[,1]

mdata<-melt(fdata)
colnames(mdata)<-c("Metabolism","Sample","Abundance")

min_abund<-1000
mdata<-mdata[which(mdata$Abundance >= min_abund),]

mdata<-merge(mdata,metadata,by="Sample",all.x=TRUE)

mdata$Zone<-"Epipelagic"
mdata$Zone[which(mdata$Depth_m > 200)] <- "Mesopelagic"
mdata$Zone[which(mdata$Depth_m > 1000)] <- "Bathypelagic"
mdata$Zone<-factor(mdata$Zone,levels=c("Epipelagic","Mesopelagic","Bathypelagic"))

summary(mdata)

metab_coloring<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69", "#A6761D","#666666")

names(metab_coloring)<-c("Amino acid metabolism","Replication and repair","Nucleotide metabolism","Metabolism of cofactors and vitamins","Signal transduction","Xenobiotics biodegradation and metabolism","Lipid metabolism","Glycan biosynthesis and metabolism","Carbohydrate metabolism","Cell growth and death","Cellular community - prokaryotes","Folding, sorting and degradation",
"Biosynthesis of other secondary metabolites","Energy metabolism","Cell motility","Metabolism of terpenoids and polyketides","Membrane transport","Metabolism of other amino acids","Transport and catabolism","Transcription","Translation")

metab_coloring<-metab_coloring[sort(as.vector(unique(mdata$Metabolism)))]

mdata$Metabolism<-factor(mdata$Metabolism,levels=sort(as.vector(unique(mdata$Metabolism))))

pdf("Profiles_Malaspina_Viral_Metabolism_Abundance_Barplot.pdf",width=16,height=16,pointsize=8)
metab_barplot<-ggplot(mdata, aes(x=Sample, y=Abundance, fill=Metabolism))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+ylab("Relative Abundance (RPKM)")+coord_flip()+theme_bw()+theme(legend.position="right",text = element_text(size = 18),legend.key.size = unit(0.5, "cm"),axis.text.x = element_text(angle = 45,hjust = 1,size=9))+scale_fill_manual(name="Metabolism",values=metab_coloring)+facet_grid(Zone ~ ., scales="free")
print(metab_barplot)
dev.off()




### 10. Metadata x Biological variables correlations (Pearson and Spearmann)
fmetadata<-subset(metadata,select=c(Temp,Conductivity,Fluo,PAR,SPAR,Turb_FTU,Sal_PSU,Salinity_WOA13,NO3_WOA13,PO4_WOA13,SiO4_WOA13,MLD,Oxygen,sigma,O2_umol_kg,O2_corr_umol_kg,O2_sat,AOU_corr_umol_kg,Chla_ugl,Fmax1_resp_prok,Fmax2_resp_euk,Fmax3_tirosina,Fmax4_triptofano,TEP,POC_uM,Turb ,pmol_leu,SE,LNA,HNA,All_BT,percentHNA,cell_size,Bacterial_cell_C,Biomass,ugC_l_d,d_1,turnover_days,HNF,low_virus,medium_virus,high_virus,all_virus,VBR))

table_files<-c("/mnt/lustre/scratch/elopez/7_normalization/host_filtered_RPKM/Filtered_by_Predicted_Host_Score-Predicted_Host_2_Phylum-RPKM.tsv","/mnt/lustre/scratch/elopez/7_normalization/virus_unfiltered_RPKM/Filtered_by_Predicted_Host_Score-family-RPKM.tsv","/mnt/lustre/scratch/elopez/7_normalization/AMG_abundance/KO_Abundance.tsv","/mnt/lustre/scratch/elopez/7_normalization/AMG_abundance/Metabolism_Abundance.tsv","/mnt/lustre/scratch/elopez/7_normalization/AMG_abundance/Pathway_Abundance.tsv")


correl_scores_matrix<-c()
for (table_file in table_files) {
	var_abund<-read.table(file=table_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
	
	rownames(var_abund)<-var_abund[,1]
	var_abund<-var_abund[,-1]
	
	for (param in colnames(fmetadata)) {
		for (variable in rownames(var_abund)) {
			var_vector<-c()
			param_vector<-c()
			for (sample_name in colnames(var_abund)) {
				sub_fmetadata<-fmetadata[which(metadata$Sample == sample_name),]
				if (nrow(sub_fmetadata) == 1) {
					param_vector<-c(param_vector,as.vector(sub_fmetadata[[param]])[1])
				} else {
					param_vector<-c(param_vector,NA)
				}
				var_vector<-c(var_vector,var_abund[variable,sample_name])
			}
		correl<-cor.test(param_vector,var_vector,method="pearson")
		scorrel<-cor.test(param_vector,var_vector,method="spearman")
		result<-c(table_file,param,variable,correl$estimate,correl$p.value,scorrel$estimate,scorrel$p.value)
		correl_scores_matrix<-rbind(correl_scores_matrix,result)
		}
	}
}

colnames(correl_scores_matrix)<-c("File","Environmental_Parameter","Metagenome_Variable","PCC","PCC_p_value","SCC","SCC_p_value")
write.table(correl_scores_matrix,file="Malaspina_Profiles_Viruses_EnvxMG_Correl_Scores.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)




### 11. Biological variables x Depth zone (Mann Whitney Tests)

# -------------- Normality tests (Shapiro-Wilk n<50) ----------------
## INPUT FILES
table_files<-c("/mnt/lustre/scratch/elopez/7_normalization/host_filtered_FPKM/Filtered_by_Predicted_Host_Score-Predicted_Host_2_Phylum-FPKM.tsv",
               "/mnt/lustre/scratch/elopez/7_normalization/virus_unfiltered_FPKM/Filtered_by_Predicted_Host_Score-family-FPKM.tsv",
               "/mnt/lustre/scratch/elopez/7_normalization/AMG_abundance/KO_Abundance.tsv",
               "/mnt/lustre/scratch/elopez/7_normalization/AMG_abundance/Metabolism_Abundance.tsv",
               "/mnt/lustre/scratch/elopez/7_normalization/AMG_abundance/Pathway_Abundance.tsv")

all_results<-as.data.frame(rbind()) # merge all input files in a data frame

## ANALYSIS
for (table_file in table_files) {
	# read file
	var_abund<-read.table(file=table_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
	
	# rename rows
	rownames(var_abund)<-var_abund[,1] # first column names to rename rows
	var_abund<-var_abund[,-1] # remove first column
	
	# modify data to make it feasible to merge them with metadata to obtain water layer
	var_abund<-as.data.frame(t(var_abund)) # traspose matrix (columns: variables; rows: samples)
	bio_vars<-as.vector(unique(colnames(var_abund)))
	var_abund$Sample<-rownames(var_abund)
	
	# merge with metadata
	var_abund<-merge(var_abund,metadata,by="Sample",all.x=TRUE)
	
	valid_zones<-as.vector(na.omit(unique(var_abund$layer)))
	
	## NORMALITY TESTS (Shapiro-Wilk n < 50; Kolmogorov-Smirnov n > 50)
	for (zone in valid_zones) {
	  for (bio_var in bio_vars) {
		vals_zone<-var_abund[which(var_abund$layer == zone),][[bio_var]] # extract Values from the Condition 
		# perform Shapiro-Wilk test (n < 50)
		if (length(vals_zone) < 50) {
		  test_result<-shapiro.test(vals_zone)
		  all_results<-rbind(all_results,c(table_file,bio_var,zone,mean(vals_zone),test_result$p.value,"Shapiro-Wilk"))
		}
	  }
	}
}

colnames(all_results)<-c("File","Variable","Group","Mean_Group","p_value","Test")

# calculate the adjusted p-value using Bonferroni method 
all_results$p_value<-as.numeric(all_results$p_value) # convert into numeric values from the p_value column


# save the results in a table
write.table(all_results,file="Malaspina_Profiles_Viruses_ZonexBiological_Vars_Normality_Test.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


# ----------------------- Fligner-Killeen -------------------------
all_results<-NULL
all_results<-as.data.frame(rbind())

for (table_file in table_files) {
	# read file
	var_abund<-read.table(file=table_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
	
	# rename rows
	rownames(var_abund)<-var_abund[,1] # first column names to rename rows
	var_abund<-var_abund[,-1] # remove first column
	
	# modify data to make it feasible to merge them with metadata to obtain water layer
	var_abund<-as.data.frame(t(var_abund)) # traspose matrix (columns: variables; rows: samples)
	bio_vars<-as.vector(unique(colnames(var_abund)))
	var_abund$Sample<-rownames(var_abund)
	
	# merge with metadata
	var_abund<-merge(var_abund,metadata,by="Sample",all.x=TRUE)
	
	valid_zones<-as.vector(na.omit(unique(var_abund$layer)))
	
	# empty vector to save pairs of values already seen
	seen_combos<-c()
	
	# for each file
	for (zoneA in valid_zones) {		
		for (zoneB in valid_zones) {
			for (bio_var in bio_vars) {
				
				if (zoneA == zoneB) { next }
				combo<-paste(zoneA,zoneB,bio_var,sep="_")
				
				if (combo %in% seen_combos) { next }
				vals_zoneA<-var_abund[which(var_abund$layer == zoneA),][[bio_var]] # extract Values from the first Condition 
				vals_zoneB<-var_abund[which(var_abund$layer == zoneB),][[bio_var]] # extract Values from the second Condition 
				
				# perform Mann-Whitney test
				test_result<-fligner.test(vals_zoneA,vals_zoneB) 
				print(test_result)
				
				# add rows with results to a dataframe (from the test, save the p-value)
				all_results<-rbind(all_results,c(table_file,bio_var,test_result$p.value)) 
				
				seen_combos<-c(seen_combos,paste(zoneB,zoneA,bio_var,sep="_")) # add to the list to not repeating the same analysis
			}
		}
	}
}

colnames(all_results)<-c("File","Variable","p_value")

# calculate the adjusted p-value using Bonferroni method 
all_results$p_value<-as.numeric(all_results$p_value) # convert into numeric values from the p_value column

# save the results in a table
write.table(all_results,file="Malaspina_Profiles_Viruses_ZonexBiological_Vars_Homoscedasticity_Test.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


# ----------------------- Mann-Whitney ----------------------------
# merge all input files in a data frame
all_results<-NULL
all_results<-as.data.frame(rbind())

for (table_file in table_files) {
	# read file
	var_abund<-read.table(file=table_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
	
	# rename rows
	rownames(var_abund)<-var_abund[,1] # first column names to rename rows
	var_abund<-var_abund[,-1] # remove first column
	
	# modify data to make it feasible to merge them with metadata to obtain water layer
	var_abund<-as.data.frame(t(var_abund)) # traspose matrix (columns: variables; rows: samples)
	bio_vars<-as.vector(unique(colnames(var_abund)))
	var_abund$Sample<-rownames(var_abund)
	
	# merge with metadata
	var_abund<-merge(var_abund,metadata,by="Sample",all.x=TRUE)
	
	valid_zones<-as.vector(na.omit(unique(var_abund$layer)))
	
	# empty vector to save pairs of values already seen
	seen_combos<-c()
	
	# for each file
	for (zoneA in valid_zones) {		
		for (zoneB in valid_zones) {
			for (bio_var in bio_vars) {
				
				if (zoneA == zoneB) { next }
				combo<-paste(zoneA,zoneB,bio_var,sep="_")
				
				if (combo %in% seen_combos) { next }
				vals_zoneA<-var_abund[which(var_abund$layer == zoneA),][[bio_var]] # extract Values from the first Condition 
				vals_zoneB<-var_abund[which(var_abund$layer == zoneB),][[bio_var]] # extract Values from the second Condition 
				
				# perform Mann-Whitney test
				test_result<-wilcox.test(vals_zoneA,vals_zoneB) 
				
				# add rows with results to a dataframe (from the test, save the p-value)
				all_results<-rbind(all_results,c(table_file,bio_var,zoneA,zoneB,mean(vals_zoneA),mean(vals_zoneB),mean(vals_zoneA)/mean(vals_zoneB),log10(mean(vals_zoneA)/mean(vals_zoneB)),test_result$p.value)) 
				
				seen_combos<-c(seen_combos,paste(zoneB,zoneA,bio_var,sep="_")) # add to the list to not repeating the same analysis
			}
		}
	}
}

colnames(all_results)<-c("File","Variable","Group_A","Group_B","Mean_A","Mean_B","Fold_Change","Log10_Fold_Change","p_value")

# calculate the adjusted p-value using Bonferroni method 
all_results$p_value<-as.numeric(all_results$p_value) # convert into numeric values from the p_value column
all_results$Adjusted_p_value<-p.adjust(all_results$p_value, method ="bonferroni")

# save the results in a table
write.table(all_results,file="Malaspina_Profiles_Viruses_ZonexBiological_Vars_Mann_Whitney.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

