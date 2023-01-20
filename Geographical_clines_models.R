# written by Juan Cubillos
#
############################################
# This script was written and compiled from two other sources
#used as reference for the present study:
#Ryan Stanley & Mallory Van Wyngaarden 
#https://github.com/rystanley/Collaborative_R_Stuff/blob/master/Mallory_HZARclines_outlier.R
#Anja Westram -Litorina hybrid zones
#https://github.com/AnjaWestram/Littorina_hybrid_zone_1/tree/master/cline%20analysis
############################################

Packages <- c("tidyverse", "tidyr", "ggplot2", "dplyr", "gridExtra","vcfR","poppr","hzar","reshape2","Hmisc","rlist","webshot2","gtExtras")
lapply(Packages, library, character.only = TRUE)

setwd("")
#Load the csv file containing the allele frequencies per population and season. Check sample file
AF<-read_csv(file.choose(),col_names = T) # Allele_Frequencies_top_20SNPs_All_populations.csv

#Add population as a Factor, here ordered greographgically from Marine to most upstream in the Hunte
AF$POP<-factor(AF$POP , levels=c("Fano","Sylt_20","Sylt_21","EMS_18","EMS_20","LG_19","VH_19","BH_18","Weser_17","Weser_20","OK_20","OK_21","HA_19","HA_20","HA_21","HA_22","HOEV","HTLS","HRM","RM","HBO","HGLA","HBW","GM_18","GM_19","HB_17","BrB","HWH","AM","FB_20", "FB_21"))
#Add Months as a factor
AF$MONTH<-factor(AF$MONTH, levels=c("February","May","June","July","October"))
#Sort the file for hzar
hzar<-AF %>% dplyr::select(CHR, POP, Distance, AF_pop,YEAR,MONTH,nSamples,TAC)

###Preparing gor HZAR
#Get the model seasonal combinations
loci<-unique(AF$CHR)# extract loci names
season<-c("May","June","July","October")
seasonal_data<-expand.grid(loci=loci,
                           season=season)
#create a list "hzar_Hunte"
hzar_Hunte<-list()
#extend the lenght of lists to the number of loci
length(hzar_Hunte)<-length(loci)
#add the loci name to each list
names(hzar_Hunte)<-loci

#create nested lists for model parameters, obs. etc into each loci
for(i in names(hzar_Hunte)){
hzar_Hunte[[i]]<-list(obs=list(),models=list(),fitRs=list(),runs=list(),analysis=list())
}
#create the lists equivalent to the seasonal sampling
for(i in names(hzar_Hunte)){
  length(hzar_Hunte[[i]]$obs)<-length(season)
}
#name the seasonal lists
for(i in names(hzar_Hunte)){
  names(hzar_Hunte[[i]]$obs)<-season
}


#create the seasonal models
#filter first to 2020 Year
hzar<-hzar %>% filter(YEAR==2020)%>% mutate(AFR=1-AF_pop) %>% dplyr::select(CHR,MONTH,POP, Distance, AF_pop, AFR, TAC,nSamples)
##Merge the seasons for parentals to a unique Month. We replace the factor February for May (all adults)
levels(hzar$MONTH)[match("February",levels(hzar$MONTH))] <- "May"

##We split the obs data between loci and seasons
obs<-split(hzar, f=list(hzar$CHR,hzar$MONTH))


#example loci G1_19692754_myosin18ab
hzar_Hunte$G1_19692754_myosin18ab$obs$May<-hzar.doMolecularData1DPops(obs$G1_19692754_myosin18ab.May$Distance,obs$G1_19692754_myosin18ab.May$AFR,obs$G1_19692754_myosin18ab.May$TAC)
hzar_Hunte$G1_19692754_myosin18ab$obs$June<-hzar.doMolecularData1DPops(obs$G1_19692754_myosin18ab.June$Distance,obs$G1_19692754_myosin18ab.June$AFR,obs$G1_19692754_myosin18ab.June$TAC)
hzar_Hunte$G1_19692754_myosin18ab$obs$July<-hzar.doMolecularData1DPops(obs$G1_19692754_myosin18ab.July$Distance,obs$G1_19692754_myosin18ab.July$AFR,obs$G1_19692754_myosin18ab.July$TAC)
hzar_Hunte$G1_19692754_myosin18ab$obs$October<-hzar.doMolecularData1DPops(obs$G1_19692754_myosin18ab.October$Distance,obs$G1_19692754_myosin18ab.October$AFR,obs$G1_19692754_myosin18ab.October$TAC)

####ALL_loci####
#"G1_21704763_atp1a1a"
hzar_Hunte$G1_21704763_atp1a1a$obs$May<-hzar.doMolecularData1DPops(obs$G1_21704763_atp1a1a.May$Distance,obs$G1_21704763_atp1a1a.May$AFR,obs$G1_21704763_atp1a1a.May$TAC)
hzar_Hunte$G1_21704763_atp1a1a$obs$June<-hzar.doMolecularData1DPops(obs$G1_21704763_atp1a1a.June$Distance,obs$G1_21704763_atp1a1a.June$AFR,obs$G1_21704763_atp1a1a.June$TAC)
hzar_Hunte$G1_21704763_atp1a1a$obs$July<-hzar.doMolecularData1DPops(obs$G1_21704763_atp1a1a.July$Distance,obs$G1_21704763_atp1a1a.July$AFR,obs$G1_21704763_atp1a1a.July$TAC)
hzar_Hunte$G1_21704763_atp1a1a$obs$October<-hzar.doMolecularData1DPops(obs$G1_21704763_atp1a1a.October$Distance,obs$G1_21704763_atp1a1a.October$AFR,obs$G1_21704763_atp1a1a.October$TAC)

#G2_405619_mucin5
hzar_Hunte$G2_405619_mucin5$obs$May<-hzar.doMolecularData1DPops(obs$G2_405619_mucin5.May$Distance,obs$G2_405619_mucin5.May$AFR,obs$G2_405619_mucin5.May$TAC)
hzar_Hunte$G2_405619_mucin5$obs$June<-hzar.doMolecularData1DPops(obs$G2_405619_mucin5.June$Distance,obs$G2_405619_mucin5.June$AFR,obs$G2_405619_mucin5.June$TAC)
hzar_Hunte$G2_405619_mucin5$obs$July<-hzar.doMolecularData1DPops(obs$G2_405619_mucin5.July$Distance,obs$G2_405619_mucin5.July$AFR,obs$G2_405619_mucin5.July$TAC)
hzar_Hunte$G2_405619_mucin5$obs$October<-hzar.doMolecularData1DPops(obs$G2_405619_mucin5.October$Distance,obs$G2_405619_mucin5.October$AFR,obs$G2_405619_mucin5.October$TAC)

#"G2_10400609_tead1db"
hzar_Hunte$G2_10400609_tead1db$obs$May<-hzar.doMolecularData1DPops(obs$G2_10400609_tead1db.May$Distance,obs$G2_10400609_tead1db.May$AFR,obs$G2_10400609_tead1db.May$TAC)
hzar_Hunte$G2_10400609_tead1db$obs$June<-hzar.doMolecularData1DPops(obs$G2_10400609_tead1db.June$Distance,obs$G2_10400609_tead1db.June$AFR,obs$G2_10400609_tead1db.June$TAC)
hzar_Hunte$G2_10400609_tead1db$obs$July<-hzar.doMolecularData1DPops(obs$G2_10400609_tead1db.July$Distance,obs$G2_10400609_tead1db.July$AFR,obs$G2_10400609_tead1db.July$TAC)
hzar_Hunte$G2_10400609_tead1db$obs$October<-hzar.doMolecularData1DPops(obs$G2_10400609_tead1db.October$Distance,obs$G2_10400609_tead1db.October$AFR,obs$G2_10400609_tead1db.October$TAC)

#G2_14538250
hzar_Hunte$G2_14538250$obs$May<-hzar.doMolecularData1DPops(obs$G2_14538250.May$Distance,obs$G2_14538250.May$AFR,obs$G2_14538250.May$TAC)
hzar_Hunte$G2_14538250$obs$June<-hzar.doMolecularData1DPops(obs$G2_14538250.June$Distance,obs$G2_14538250.June$AFR,obs$G2_14538250.June$TAC)
hzar_Hunte$G2_14538250$obs$July<-hzar.doMolecularData1DPops(obs$G2_14538250.July$Distance,obs$G2_14538250.July$AFR,obs$G2_14538250.July$TAC)
hzar_Hunte$G2_14538250$obs$October<-hzar.doMolecularData1DPops(obs$G2_14538250.October$Distance,obs$G2_14538250.October$AFR,obs$G2_14538250.October$TAC)

#"G4_12806385_EDA"
hzar_Hunte$G4_12806385_EDA$obs$May<-hzar.doMolecularData1DPops(obs$G4_12806385_EDA.May$Distance,obs$G4_12806385_EDA.May$AFR,obs$G4_12806385_EDA.May$TAC)
hzar_Hunte$G4_12806385_EDA$obs$June<-hzar.doMolecularData1DPops(obs$G4_12806385_EDA.June$Distance,obs$G4_12806385_EDA.June$AFR,obs$G4_12806385_EDA.June$TAC)
hzar_Hunte$G4_12806385_EDA$obs$July<-hzar.doMolecularData1DPops(obs$G4_12806385_EDA.July$Distance,obs$G4_12806385_EDA.July$AFR,obs$G4_12806385_EDA.July$TAC)
hzar_Hunte$G4_12806385_EDA$obs$October<-hzar.doMolecularData1DPops(obs$G4_12806385_EDA.October$Distance,obs$G4_12806385_EDA.October$AFR,obs$G4_12806385_EDA.October$TAC)

#"G4_19879931" 
hzar_Hunte$G4_19879931$obs$May<-hzar.doMolecularData1DPops(obs$G4_19879931.May$Distance,obs$G4_19879931.May$AFR,obs$G4_19879931.May$TAC)
hzar_Hunte$G4_19879931$obs$June<-hzar.doMolecularData1DPops(obs$G4_19879931.June$Distance,obs$G4_19879931.June$AFR,obs$G4_19879931.June$TAC)
hzar_Hunte$G4_19879931$obs$July<-hzar.doMolecularData1DPops(obs$G4_19879931.July$Distance,obs$G4_19879931.July$AFR,obs$G4_19879931.July$TAC)
hzar_Hunte$G4_19879931$obs$October<-hzar.doMolecularData1DPops(obs$G4_19879931.October$Distance,obs$G4_19879931.October$AFR,obs$G4_19879931.October$TAC)

#"G7_1829245"   
hzar_Hunte$G7_1829245$obs$May<-hzar.doMolecularData1DPops(obs$G7_1829245.May$Distance,obs$G7_1829245.May$AFR,obs$G7_1829245.May$TAC)
hzar_Hunte$G7_1829245$obs$June<-hzar.doMolecularData1DPops(obs$G7_1829245.June$Distance,obs$G7_1829245.June$AFR,obs$G7_1829245.June$TAC)
hzar_Hunte$G7_1829245$obs$July<-hzar.doMolecularData1DPops(obs$G7_1829245.July$Distance,obs$G7_1829245.July$AFR,obs$G7_1829245.July$TAC)
hzar_Hunte$G7_1829245$obs$October<-hzar.doMolecularData1DPops(obs$G7_1829245.October$Distance,obs$G7_1829245.October$AFR,obs$G7_1829245.October$TAC)

#"G7_20827694_mtnr1bb"
hzar_Hunte$G7_20827694_mtnr1bb$obs$May<-hzar.doMolecularData1DPops(obs$G7_20827694_mtnr1bb.May$Distance,obs$G7_20827694_mtnr1bb.May$AFR,obs$G7_20827694_mtnr1bb.May$TAC)
hzar_Hunte$G7_20827694_mtnr1bb$obs$June<-hzar.doMolecularData1DPops(obs$G7_20827694_mtnr1bb.June$Distance,obs$G7_20827694_mtnr1bb.June$AFR,obs$G7_20827694_mtnr1bb.June$TAC)
hzar_Hunte$G7_20827694_mtnr1bb$obs$July<-hzar.doMolecularData1DPops(obs$G7_20827694_mtnr1bb.July$Distance,obs$G7_20827694_mtnr1bb.July$AFR,obs$G7_20827694_mtnr1bb.July$TAC)
hzar_Hunte$G7_20827694_mtnr1bb$obs$October<-hzar.doMolecularData1DPops(obs$G7_20827694_mtnr1bb.October$Distance,obs$G7_20827694_mtnr1bb.October$AFR,obs$G7_20827694_mtnr1bb.October$TAC)

#"G8_1302871" 
hzar_Hunte$G8_1302871$obs$May<-hzar.doMolecularData1DPops(obs$G8_1302871.May$Distance,obs$G8_1302871.May$AFR,obs$G8_1302871.May$TAC)
hzar_Hunte$G8_1302871$obs$June<-hzar.doMolecularData1DPops(obs$G8_1302871.June$Distance,obs$G8_1302871.June$AFR,obs$G8_1302871.June$TAC)
hzar_Hunte$G8_1302871$obs$July<-hzar.doMolecularData1DPops(obs$G8_1302871.July$Distance,obs$G8_1302871.July$AFR,obs$G8_1302871.July$TAC)
hzar_Hunte$G8_1302871$obs$October<-hzar.doMolecularData1DPops(obs$G8_1302871.October$Distance,obs$G8_1302871.October$AFR,obs$G8_1302871.October$TAC)

#"G9_13197858" 
hzar_Hunte$G9_13197858$obs$May<-hzar.doMolecularData1DPops(obs$G9_13197858.May$Distance,obs$G9_13197858.May$AFR,obs$G9_13197858.May$TAC)
hzar_Hunte$G9_13197858$obs$June<-hzar.doMolecularData1DPops(obs$G9_13197858.June$Distance,obs$G9_13197858.June$AFR,obs$G9_13197858.June$TAC)
hzar_Hunte$G9_13197858$obs$July<-hzar.doMolecularData1DPops(obs$G9_13197858.July$Distance,obs$G9_13197858.July$AFR,obs$G9_13197858.July$TAC)
hzar_Hunte$G9_13197858$obs$October<-hzar.doMolecularData1DPops(obs$G9_13197858.October$Distance,obs$G9_13197858.October$AFR,obs$G9_13197858.October$TAC)

#"G10_4238059"   
hzar_Hunte$G10_4238059$obs$May<-hzar.doMolecularData1DPops(obs$G10_4238059.May$Distance,obs$G10_4238059.May$AFR,obs$G10_4238059.May$TAC)
hzar_Hunte$G10_4238059$obs$June<-hzar.doMolecularData1DPops(obs$G10_4238059.June$Distance,obs$G10_4238059.June$AFR,obs$G10_4238059.June$TAC)
hzar_Hunte$G10_4238059$obs$July<-hzar.doMolecularData1DPops(obs$G10_4238059.July$Distance,obs$G10_4238059.July$AFR,obs$G10_4238059.July$TAC)
hzar_Hunte$G10_4238059$obs$October<-hzar.doMolecularData1DPops(obs$G10_4238059.October$Distance,obs$G10_4238059.October$AFR,obs$G10_4238059.October$TAC)

#"G12_11155054_arhgef25"
hzar_Hunte$G12_11155054_arhgef25$obs$May<-hzar.doMolecularData1DPops(obs$G12_11155054_arhgef25.May$Distance,obs$G12_11155054_arhgef25.May$AFR,obs$G12_11155054_arhgef25.May$TAC)
hzar_Hunte$G12_11155054_arhgef25$obs$June<-hzar.doMolecularData1DPops(obs$G12_11155054_arhgef25.June$Distance,obs$G12_11155054_arhgef25.June$AFR,obs$G12_11155054_arhgef25.June$TAC)
hzar_Hunte$G12_11155054_arhgef25$obs$July<-hzar.doMolecularData1DPops(obs$G12_11155054_arhgef25.July$Distance,obs$G12_11155054_arhgef25.July$AFR,obs$G12_11155054_arhgef25.July$TAC)
hzar_Hunte$G12_11155054_arhgef25$obs$October<-hzar.doMolecularData1DPops(obs$G12_11155054_arhgef25.October$Distance,obs$G12_11155054_arhgef25.October$AFR,obs$G12_11155054_arhgef25.October$TAC)

#"G13_8453167_EPX"
hzar_Hunte$G13_8453167_EPX$obs$May<-hzar.doMolecularData1DPops(obs$G13_8453167_EPX.May$Distance,obs$G13_8453167_EPX.May$AFR,obs$G13_8453167_EPX.May$TAC)
hzar_Hunte$G13_8453167_EPX$obs$June<-hzar.doMolecularData1DPops(obs$G13_8453167_EPX.June$Distance,obs$G13_8453167_EPX.June$AFR,obs$G13_8453167_EPX.June$TAC)
hzar_Hunte$G13_8453167_EPX$obs$July<-hzar.doMolecularData1DPops(obs$G13_8453167_EPX.July$Distance,obs$G13_8453167_EPX.July$AFR,obs$G13_8453167_EPX.July$TAC)
hzar_Hunte$G13_8453167_EPX$obs$October<-hzar.doMolecularData1DPops(obs$G13_8453167_EPX.October$Distance,obs$G13_8453167_EPX.October$AFR,obs$G13_8453167_EPX.October$TAC)

#"G14_15149501_SPTLC1"
hzar_Hunte$G14_15149501_SPTLC1$obs$May<-hzar.doMolecularData1DPops(obs$G14_15149501_SPTLC1.May$Distance,obs$G14_15149501_SPTLC1.May$AFR,obs$G14_15149501_SPTLC1.May$TAC)
hzar_Hunte$G14_15149501_SPTLC1$obs$June<-hzar.doMolecularData1DPops(obs$G14_15149501_SPTLC1.June$Distance,obs$G14_15149501_SPTLC1.June$AFR,obs$G14_15149501_SPTLC1.June$TAC)
hzar_Hunte$G14_15149501_SPTLC1$obs$July<-hzar.doMolecularData1DPops(obs$G14_15149501_SPTLC1.July$Distance,obs$G14_15149501_SPTLC1.July$AFR,obs$G14_15149501_SPTLC1.July$TAC)
hzar_Hunte$G14_15149501_SPTLC1$obs$October<-hzar.doMolecularData1DPops(obs$G14_15149501_SPTLC1.October$Distance,obs$G14_15149501_SPTLC1.October$AFR,obs$G14_15149501_SPTLC1.October$TAC)

#"G15_7634626"                
hzar_Hunte$G15_7634626$obs$May<-hzar.doMolecularData1DPops(obs$G15_7634626.May$Distance,obs$G15_7634626.May$AFR,obs$G15_7634626.May$TAC)
hzar_Hunte$G15_7634626$obs$June<-hzar.doMolecularData1DPops(obs$G15_7634626.June$Distance,obs$G15_7634626.June$AFR,obs$G15_7634626.June$TAC)
hzar_Hunte$G15_7634626$obs$July<-hzar.doMolecularData1DPops(obs$G15_7634626.July$Distance,obs$G15_7634626.July$AFR,obs$G15_7634626.July$TAC)
hzar_Hunte$G15_7634626$obs$October<-hzar.doMolecularData1DPops(obs$G15_7634626.October$Distance,obs$G15_7634626.October$AFR,obs$G15_7634626.October$TAC)

#"G17_6663596_close_to_kif21b" 
hzar_Hunte$G17_6663596_close_to_kif21b$obs$May<-hzar.doMolecularData1DPops(obs$G17_6663596_close_to_kif21b.May$Distance,obs$G17_6663596_close_to_kif21b.May$AFR,obs$G17_6663596_close_to_kif21b.May$TAC)
hzar_Hunte$G17_6663596_close_to_kif21b$obs$June<-hzar.doMolecularData1DPops(obs$G17_6663596_close_to_kif21b.June$Distance,obs$G17_6663596_close_to_kif21b.June$AFR,obs$G17_6663596_close_to_kif21b.June$TAC)
hzar_Hunte$G17_6663596_close_to_kif21b$obs$July<-hzar.doMolecularData1DPops(obs$G17_6663596_close_to_kif21b.July$Distance,obs$G17_6663596_close_to_kif21b.July$AFR,obs$G17_6663596_close_to_kif21b.July$TAC)
hzar_Hunte$G17_6663596_close_to_kif21b$obs$October<-hzar.doMolecularData1DPops(obs$G17_6663596_close_to_kif21b.October$Distance,obs$G17_6663596_close_to_kif21b.October$AFR,obs$G17_6663596_close_to_kif21b.October$TAC)

#"G18_10525993_zc3h14"        
hzar_Hunte$G18_10525993_zc3h14$obs$May<-hzar.doMolecularData1DPops(obs$G18_10525993_zc3h14.May$Distance,obs$G18_10525993_zc3h14.May$AFR,obs$G18_10525993_zc3h14.May$TAC)
hzar_Hunte$G18_10525993_zc3h14$obs$June<-hzar.doMolecularData1DPops(obs$G18_10525993_zc3h14.June$Distance,obs$G18_10525993_zc3h14.June$AFR,obs$G18_10525993_zc3h14.June$TAC)
hzar_Hunte$G18_10525993_zc3h14$obs$July<-hzar.doMolecularData1DPops(obs$G18_10525993_zc3h14.July$Distance,obs$G18_10525993_zc3h14.July$AFR,obs$G18_10525993_zc3h14.July$TAC)
hzar_Hunte$G18_10525993_zc3h14$obs$October<-hzar.doMolecularData1DPops(obs$G18_10525993_zc3h14.October$Distance,obs$G18_10525993_zc3h14.October$AFR,obs$G18_10525993_zc3h14.October$TAC)

#"G20_13028350"
hzar_Hunte$G20_13028350$obs$May<-hzar.doMolecularData1DPops(obs$G20_13028350.May$Distance,obs$G20_13028350.May$AFR,obs$G20_13028350.May$TAC)
hzar_Hunte$G20_13028350$obs$June<-hzar.doMolecularData1DPops(obs$G20_13028350.June$Distance,obs$G20_13028350.June$AFR,obs$G20_13028350.June$TAC)
hzar_Hunte$G20_13028350$obs$July<-hzar.doMolecularData1DPops(obs$G20_13028350.July$Distance,obs$G20_13028350.July$AFR,obs$G20_13028350.July$TAC)
hzar_Hunte$G20_13028350$obs$October<-hzar.doMolecularData1DPops(obs$G20_13028350.October$Distance,obs$G20_13028350.October$AFR,obs$G20_13028350.October$TAC)


####MODELS#####
model<-c(paste0("model ",1:15))
s<-c("none","fixed","free")
t<-c("none", "right", "left", "mirror", "both")
models_hzar<-expand.grid(scaling=c("none","fixed","free"),
                         tails=c("none","left","right","mirror","both"))
models_hzar<-cbind(model, models_hzar)
##to export a formated table for models
library(gtExtras)

models_table_summary<-models_hzar %>% 
  gt() %>% 
  gt_theme_nytimes() %>% 
  cols_width(
    model ~ px(120),
    scaling ~ px(100),
    tails ~ px(100)) %>% 
  tab_header(title = "Cline Models fitted with variable parameters")

gtsave(models_table_summary, "Permutations/models_table_summary.png")

#create the lists equivalent to the seasonal models
for(i in names(hzar_Hunte)){
  length(hzar_Hunte[[i]]$models)<-length(season)
}
#name the seasonal lists in models
for(i in names(hzar_Hunte)){
  names(hzar_Hunte[[i]]$models)<-season
}

#create the number of models with different scaling and tails argument options 
for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$models)){
  hzar_Hunte[[i]]$models[[j]]<-list("model_1"=list(),"model_2"=list(),"model_3"=list(),"model_4"=list(),"model_5"=list(),"model_6"=list(), "model_7"=list(),
                                 "model_8"=list(),"model_9"=list(),"model_10"=list(),"model_11"=list(),"model_12"=list(),"model_13"=list(),"model_14"=list(),"model_15"=list())
  }
}

####MODELS commands####
#May
#model_1
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_1 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="none", tails="none")
}
#model_2
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_2 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="fixed", tails="none")
}
#model_3
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_3 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="free", tails="none")
}
#model_4
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_4 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="none", tails="left")
}
#model_5
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_5 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="fixed", tails="left")
}
#model_6
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_6 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="free", tails="left")
}
#model_7
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_7 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="none", tails="right")
}
#model_8
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_8 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="fixed", tails="right")
}
#model_9
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_9 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="free", tails="right")
}
#model_10
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_10 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="none", tails="mirror")
}
#model_11
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_11 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="fixed", tails="mirror")
}
#model_12
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_12 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="free", tails="mirror")
}
#model_13
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_13 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="none", tails="both")
}
#model_14
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_14 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="fixed", tails="both")
}
#model_15
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$May$model_15 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                             scaling="free", tails="both")
}


##June
#model_1
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_1 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="none", tails="none")
}
#model_2
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_2 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="fixed", tails="none")
}
#model_3
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_3 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="free", tails="none")
}
#model_4
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_4 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="none", tails="left")
}
#model_5
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_5 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="fixed", tails="left")
}
#model_6
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_6 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="free", tails="left")
}
#model_7
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_7 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="none", tails="right")
}
#model_8
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_8 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="fixed", tails="right")
}
#model_9
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_9 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                             scaling="free", tails="right")
}
#model_10
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_10 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                              scaling="none", tails="mirror")
}
#model_11
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_11 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                              scaling="fixed", tails="mirror")
}
#model_12
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_12 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                              scaling="free", tails="mirror")
}
#model_13
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_13 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                              scaling="none", tails="both")
}
#model_14
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_14 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                              scaling="fixed", tails="both")
}
#model_15
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$June$model_15 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$June, 
                                                              scaling="free", tails="both")
}

##July
#model_1
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_1 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="none", tails="none")
}
#model_2
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_2 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="fixed", tails="none")
}
#model_3
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_3 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="free", tails="none")
}
#model_4
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_4 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="none", tails="left")
}
#model_5
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_5 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="fixed", tails="left")
}
#model_6
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_6 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="free", tails="left")
}
#model_7
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_7 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="none", tails="right")
}
#model_8
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_8 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="fixed", tails="right")
}
#model_9
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_9 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                             scaling="free", tails="right")
}
#model_10
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_10 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                              scaling="none", tails="mirror")
}
#model_11
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_11 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                              scaling="fixed", tails="mirror")
}
#model_12
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_12 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                              scaling="free", tails="mirror")
}
#model_13
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_13 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                              scaling="none", tails="both")
}
#model_14
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_14 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                              scaling="fixed", tails="both")
}
#model_15
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$July$model_15 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$July, 
                                                              scaling="free", tails="both")
}

#October
#model_1
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_1 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="none", tails="none")
}
#model_2
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_2 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="fixed", tails="none")
}
#model_3
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_3 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="free", tails="none")
}
#model_4
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_4 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="none", tails="left")
}
#model_5
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_5 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="fixed", tails="left")
}
#model_6
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_6 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="free", tails="left")
}
#model_7
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_7 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="none", tails="right")
}
#model_8
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_8 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="fixed", tails="right")
}
#model_9
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_9 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                             scaling="free", tails="right")
}
#model_10
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_10 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                              scaling="none", tails="mirror")
}
#model_11
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_11 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                              scaling="fixed", tails="mirror")
}
#model_12
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_12 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                              scaling="free", tails="mirror")
}
#model_13
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_13 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                              scaling="none", tails="both")
}
#model_14
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_14 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$October, 
                                                              scaling="fixed", tails="both")
}
#model_15
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$models$October$model_15 <- hzar.makeCline1DFreq(data=hzar_Hunte[[i]]$obs$May, 
                                                              scaling="free", tails="both")
}




#####

#Check the model parameters
print(hzar_Hunte$G1_19692754_myosin18ab$models)

#modify all models to focus on observed region 
#data collected between 0 and 4400km
for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$models)){
  hzar_Hunte[[i]]$models[[j]] <- sapply(hzar_Hunte[[i]]$models[[j]], hzar.model.addBoxReq,
                                 0, 300, simplify=FALSE)
  }}

#Check the updated model parameters for lower= 0 and upper=300
print(hzar_Hunte$G1_19692754_myosin18ab$models)

##create init and chains lists inside fitRs
for(i in names(hzar_Hunte)){
  hzar_Hunte[[i]]$fitRs<-list("init"=list(),"chains"=list())
  }

#create the lists equivalent to the seasonal models withing init (parameters) and chains (mcmc runs params)
for(i in names(hzar_Hunte)){
  length(hzar_Hunte[[i]]$fitRs$init)<-length(season)
  length(hzar_Hunte[[i]]$fitRs$chains)<-length(season)
}
#name the seasonal lists in models with the seasonals
for(i in names(hzar_Hunte)){
  names(hzar_Hunte[[i]]$fitRs$init)<-season
  names(hzar_Hunte[[i]]$fitRs$chains)<-season
}

###compile models to prepare for fitting --> creates hzar.fitRequest from each 
#clineMetalModel object

for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$models)){
    hzar_Hunte[[i]]$fitRs$init[[j]] <- sapply(hzar_Hunte[[i]]$models[[j]], 
                                     hzar.first.fitRequest.old.ML,
                                     obsData = hzar_Hunte[[i]]$obs[[j]],
                                     verbose=FALSE,
                                     simplify=FALSE)


}}

###Running the models####

## A typical chain length.  The default value is 1e5 However in previous analysis we found convergence around 2e5 and 6e5 for exceptional cases).
chainLength=2e5;

## Make each model run off a separate seed
#random seeds can be calculated as floor(runif(24, min=0, max=1000))
mainSeed=
  list(A=c(765,41,38,435,885,229),
       B=c(400,552,69,210,922,299),
       C=c(602,551,689,688,46,353),
       D=c(973,411,132,146,594,972))

###update models chain length and burnin
#update settings for the fitter using chainLength and mainSeed created before
for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$fitRs$init)){
    for (k in names(hzar_Hunte[[i]]$fitRs$init[[j]])){
    hzar_Hunte[[i]]$fitRs$init[[j]][[k]]$mcmcParam$chainLength <- chainLength}}}

for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$fitRs$init)){
    for (k in names(hzar_Hunte[[i]]$fitRs$init[[j]])){
    hzar_Hunte[[i]]$fitRs$init[[j]][[k]]$mcmcParam$burnin <- chainLength %/% 10}}}

#check fit request settings 
print(hzar_Hunte$G1_19692754_myosin18ab$fitRs$init$May$model_5)

#replicate each fit request 3 times, keeping the original seeds 
#while switching to a new seed channel
#45 total fit requests - 15 models, 3 times each
for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$fitRs$init)){
    for (k in names(hzar_Hunte[[i]]$fitRs$init[[j]])){
  hzar_Hunte[[i]]$fitRs$chains[[j]][[k]] <- hzar.multiFitRequest(hzar_Hunte[[i]]$fitRs$init[[j]][[k]],
                                                     each=3,
                                                     baseSeed=NULL)
}}}

#create the lists equivalent to the seasonal models withing init (parameters) and chains (mcmc runs params)
for(i in names(hzar_Hunte)){
  length(hzar_Hunte[[i]]$runs)<-length(season)
}
#name the seasonal lists in models with the seasonals
for(i in names(hzar_Hunte)){
  names(hzar_Hunte[[i]]$runs)<-season
}

for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$runs)){
    hzar_Hunte[[i]]$runs[[j]]<-list("model_1"=list(),"model_2"=list(),"model_3"=list(),"model_4"=list(),"model_5"=list(),"model_6"=list(), "model_7"=list(),
                                      "model_8"=list(),"model_9"=list(),"model_10"=list(),"model_11"=list(),"model_12"=list(),"model_13"=list(),"model_14"=list(),"model_15"=list())
  }
}


##have 3420 fit requests - models 15, each with 3 chain, 
#running the chain 3 times - 3420 total runs? - THIS WILL TAKE AROUND 24 HOURS!!
for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$fitRs$init)){
    for (k in names(hzar_Hunte[[i]]$fitRs$init[[j]])){
      hzar_Hunte[[i]]$runs[[j]][[k]]$doSeq <- lapply(hzar_Hunte[[i]]$fitRs$chains[[j]][[k]], hzar.chain.doSeq,
                                                       count=3)
    }}}


#save the results and the hzar_Hunte lists into an RDS file.
#list.save(hzar_Hunte, "hzar_Hunte_list.rds")
#restore the file with the following command if needed
#readRDS(, refhook = NULL)

##CHeck for model convergence
#this command will display all the runs parameters and results
hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_1$doSeq[1:3]
hzar_Hunte_2$G1_21704763_atp1a1a$runs$May$model_1$doSeq[1:3]

#this command summarizes all the 9 runs (3 runs each with 3 LL)
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_1$doSeq[1:3],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#summary(do.call(mcmc.list, lapply(hzar_Hunte_2$G1_21704763_atp1a1a$runs$May$model_1$doSeq[1:3],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))

summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_2$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_3$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_4$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_5$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_6$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_7$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_8$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_9$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_10$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_11$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_12$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_13$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_14$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))
summary(do.call(mcmc.list, lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_15$doSeq[1:3],
                                  function(x) plot(hzar.mcmc.bindLL(x[[3]])))))


###ANALYSIS####
#start aggregation of data for analysis
#CREATE SEASONAL LISTS for the analysis lists inside each loci

hzar_Hunte<-readRDS("/home/juan/Documents/Sticklebacks/GTS_Pool1/ADMIXTURE/plink_linux_x86_64_20210606/Course/Allele_Freq/All_Runs/Allele_frequecies_per_Season_Pop/hzar_Hunte_list.rds")
hzar_Hunte_2<-readRDS("/home/juan/Documents/Sticklebacks/GTS_Pool1/ADMIXTURE/plink_linux_x86_64_20210606/Course/Allele_Freq/All_Runs/Allele_frequecies_per_Season_Pop/hzar_Hunte_list.rds")

#create the lists equivalent to the seasonal models within init (parameters) and chains (mcmc runs params)
for(i in names(hzar_Hunte)){
  length(hzar_Hunte[[i]]$analysis)<-length(season)
}
#name the seasonal lists in models with the seasonal name
for(i in names(hzar_Hunte)){
  names(hzar_Hunte[[i]]$analysis)<-season
}

##now add all the other models where the multiple fits will be grouped
for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$analysis)){
    hzar_Hunte[[i]]$analysis[[j]]$initDGs<-list("nullModel"=list(),"model_1"=list(),"model_2"=list(),"model_3"=list(),"model_4"=list(),"model_5"=list(),"model_6"=list(), "model_7"=list(),
                                    "model_8"=list(),"model_9"=list(),"model_10"=list(),"model_11"=list(),"model_12"=list(),"model_13"=list(),"model_14"=list(),"model_15"=list())
  }
}

##create a data group for the null model
for(i in names(hzar_Hunte)){
  for (j in names(hzar_Hunte[[i]]$analysis)){
    hzar_Hunte[[i]]$analysis[[j]]$initDGs$nullModel <- hzar.dataGroup.null(hzar_Hunte[[i]]$obs[[j]])
  }}


#create a model data group for each model from the initial runs
for(i in names(hzar_Hunte)){ # for the loci name
  for (j in names(hzar_Hunte[[i]]$runs)){ #Season name
    for (k in names(hzar_Hunte[[i]]$runs[[j]])){ # model number
      for(l in names(hzar_Hunte[[i]]$runs[[j]][[k]]$doSeq)){ #for the run
       hzar_Hunte[[i]]$analysis[[j]]$initDGs[[k]] <- hzar.dataGroup.add(lapply(hzar_Hunte[[i]]$runs[[j]][[k]][[doSeq]],
                                                                             hzar.dataGroup.add),
                                                                      hzar_Hunte[[i]]$analysis[[j]]$initDGs[[k]])
     }
    }
  }
}


for(i in names(hzar_Hunte)){ # for the loci name
  for (j in names(hzar_Hunte[[i]]$runs)){ #Season name
    for (k in names(hzar_Hunte[[i]]$runs[[j]])){ # model number
      for(l in names(hzar_Hunte[[i]]$runs[[j]][[k]]$doSeq)){ #for the run
        # for(m in names(hzar_Hunte[[i]]$runs[[j]][[k]]$doSeq[[l]])){
        hzar_Hunte[[i]]$analysis[[j]]$initDGs[[k]] <- hzar.dataGroup.add(hzar_Hunte[[i]]$runs[[j]][[k]][[doSeq]])
      }
    }
  }
}

####create a model data group for each model from the initial runs####
#G1_19692754_myosin18ab####
#May#
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_19692754_myosin18ab$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_19692754_myosin18ab$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#G1_21704763_atp1a1a####
#May#
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#June#
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G1_21704763_atp1a1a$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G1_21704763_atp1a1a$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G2_405619_mucin5####
#May#
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_405619_mucin5$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_405619_mucin5$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#G2_10400609_tead1db####
#May#
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_10400609_tead1db$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_10400609_tead1db$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#G2_14538250####
#May#
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G2_14538250$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G2_14538250$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#G4_12806385_EDA####
#May#
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_12806385_EDA$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_12806385_EDA$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G4_19879931####
#May#
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G4_19879931$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G4_19879931$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G7_1829245####
#May#
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_1829245$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_1829245$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G7_20827694_mtnr1bb####
#May#
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G7_20827694_mtnr1bb$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G7_20827694_mtnr1bb$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G8_1302871####
#May#
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G8_1302871$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G8_1302871$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G9_13197858####
#May#
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G9_13197858$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G9_13197858$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)



#G10_4238059####
#May#
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G10_4238059$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G10_4238059$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G12_11155054_arhgef25####
#May#
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G12_11155054_arhgef25$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G12_11155054_arhgef25$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G13_8453167_EPX####
#May#
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G13_8453167_EPX$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G13_8453167_EPX$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G14_15149501_SPTLC1####
#May#
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G14_15149501_SPTLC1$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G14_15149501_SPTLC1$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G15_7634626####
#May#
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G15_7634626$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G15_7634626$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G17_6663596_close_to_kif21b####
#May#
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G17_6663596_close_to_kif21b$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G17_6663596_close_to_kif21b$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)


#G18_10525993_zc3h14####
#May#
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G18_10525993_zc3h14$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G18_10525993_zc3h14$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)



#G20_13028350####
#May#
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$May$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$May$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#June#
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$June$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$June$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#July#
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$July$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$July$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#October#
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_1 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_1$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_2 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_2$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_3 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_3$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_4 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_4$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_5 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_5$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_6 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_6$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_7 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_7$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_8 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_8$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_9 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_9$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_10 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_10$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_11 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_11$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_12 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_12$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_13 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_13$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_14 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_14$doSeq[[1]], fitRequestL = list(), doPar = TRUE)
hzar_Hunte$G20_13028350$analysis$October$initDGs$model_15 <- hzar.dataGroup.add(hzar_Hunte$G20_13028350$runs$October$model_15$doSeq[[1]], fitRequestL = list(), doPar = TRUE)

#####

#list.save(hzar_Hunte, "hzar_Hunte_list.rds")

##create a hzar.obsDataGroup object from the 15 hzar.dataGroup just created, copying the naming scheme
hzar_Hunte$G1_19692754_myosin18ab$analysis$May$oDG <- hzar.make.obsDataGroup(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs)

hzar_Hunte$G1_19692754_myosin18ab$analysis$May$oDG <- hzar.copyModelLabels(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$initDGs,
                                                                           hzar_Hunte$G1_19692754_myosin18ab$analysis$May$oDG)

##create a hzar.obsDataGroup object from the 15 hzar.dataGroup just created, copying the naming scheme
for(i in names(hzar_Hunte)){# for the loci name
  for (j in names(hzar_Hunte[[i]]$analysis)){ #Season name
    hzar_Hunte[[i]]$analysis[[j]]$oDG <- hzar.make.obsDataGroup(lapply(hzar_Hunte[[i]]$analysis[[j]]$initDGs,
                                                                       hzar.dataGroup.add),
                                                                hzar_Hunte[[i]]$analysis[[j]]$oDG)
    }}

for(i in names(hzar_Hunte)){# for the loci name
  for (j in names(hzar_Hunte[[i]]$analysis)){ #Season name
    hzar_Hunte[[i]]$analysis[[j]]$oDG <- hzar.copyModelLabels(hzar_Hunte[[i]]$analysis[[j]]$initDGs,
                                                                     hzar_Hunte[[i]]$analysis[[j]]$oDG)
  }
}


#list.save(hzar_Hunte, "hzar_Hunte_list.rds")


##convert all 36 runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object
for(i in names(hzar_Hunte)){# for the loci name
  for (j in names(hzar_Hunte[[i]]$runs)){ #Season name
    for (k in names(hzar_Hunte[[i]]$runs[[j]])){# model number
      hzar_Hunte[[i]]$analysis[[j]]$oDG <- hzar.make.obsDataGroup(lapply(hzar_Hunte[[i]]$runs[[j]][[k]]$doSeq,
                                                                           hzar.dataGroup.add),
                                                                    hzar_Hunte[[i]]$analysis[[j]]$oDG)
      }}}

#list.save(hzar_Hunte, "hzar_Hunte_list.rds")

#for(i in names(hzar_Hunte$G1_19692754_myosin18ab$runs$May)){
 # hzar_Hunte$G1_19692754_myosin18ab$analysis$May$oDG <- hzar.make.obsDataGroup(lapply(hzar_Hunte$G1_19692754_myosin18ab$runs$May[[i]]$doSeq, 
 #                                                             hzar.dataGroup.add),
  #                                                            hzar_Hunte$G1_19692754_myosin18ab$analysis$May$oDG)
#}

#check to make sure there are only 5 hzar.dataGroup objects
print(summary(hzar_Hunte$G18_10525993_zc3h14$analysis$May$oDG$data.groups))

#compare the 15 cline models to the null model graphically
hzar.plot.cline(hzar_Hunte$G1_19692754_myosin18ab$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G1_21704763_atp1a1a$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G2_405619_mucin5$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G2_10400609_tead1db$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G2_14538250$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G4_12806385_EDA$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G4_19879931$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G7_1829245$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G7_20827694_mtnr1bb$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G8_1302871$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G9_13197858$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G10_4238059$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G12_11155054_arhgef25$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G13_8453167_EPX$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G14_15149501_SPTLC1$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G15_7634626$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G17_6663596_close_to_kif21b$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G18_10525993_zc3h14$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G20_13028350$analysis$June$oDG)
hzar.plot.cline(hzar_Hunte$G20_13028350$analysis$July$oDG)
hzar.plot.cline(hzar_Hunte$G20_13028350$analysis$October$oDG)


#model selection based on AICc scores
#print the model with the minimum AICc score
for(i in names(hzar_Hunte)){# for the loci name
  for (j in names(hzar_Hunte[[i]]$runs)){ #Season name
    for (k in names(hzar_Hunte[[i]]$runs[[j]])){# model number
      hzar_Hunte[[i]]$analysis[[j]]$AICcTable <-hzar.AICc.hzar.obsDataGroup(hzar_Hunte[[i]]$analysis[[j]]$oDG)
    }}}

#print(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$AICcTable <- 
#        hzar.AICc.hzar.obsDataGroup(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$oDG))

#list.save(hzar_Hunte, "hzar_Hunte_list.rds")

#print the model with the minimum AICc score
for(i in names(hzar_Hunte)){# for the loci name
  for (j in names(hzar_Hunte[[i]]$runs)){ #Season name
    for (k in names(hzar_Hunte[[i]]$runs[[j]])){# model number
      print(hzar_Hunte[[i]]$analysis[[j]]$model.name<-rownames(hzar_Hunte[[i]]$analysis[[j]]$AICcTable)
            [[which.min(hzar_Hunte[[i]]$analysis[[j]]$AICcTable$AICc)]])
      }}}

#print(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.name <-
#        rownames(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$AICcTable)[[which.min(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$AICcTable$AICc )]])
#[1] "model3"

#Extract the hzar.dataGroup object for the selected model
for(i in names(hzar_Hunte)){# for the loci name
  for (j in names(hzar_Hunte[[i]]$analysis)){ #Season name
    hzar_Hunte[[i]]$analysis[[j]]$model.selected<-hzar_Hunte[[i]]$analysis[[j]]$oDG$data.groups[[hzar_Hunte[[i]]$analysis[[j]]$model.name]]
            }}


#hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.selected <-
#  hzar_Hunte$G1_19692754_myosin18ab$analysis$May$oDG$data.groups[[hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.name]]

#look at the variation in parameters for the selected model
print(hzar.getLLCutParam(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.selected,
                         names(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.selected$data.param)))

####Print Params####

#print the cline width for the selected model
for(i in names(hzar_Hunte)){# for the loci name
  for (j in names(hzar_Hunte[[i]]$analysis)){ #Season name
    hzar_Hunte[[i]]$analysis[[j]]$modeldetails<-print(hzar.get.ML.cline(hzar_Hunte[[i]]$analysis[[j]]$model.selected))
  }}


hzar_Hunte$G1_19692754_myosin18ab$analysis$May$modeldetails <- print(hzar.get.ML.cline(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.selected))
print(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$modeldetails$param.all$width)

print(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$modeldetails$logLike)

#list.save(hzar_Hunte, "hzar_Hunte_list.rds")

#plot the maximum likelihood cline for the selected model
G1_19_May<-hzar.plot.cline(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.selected)

#plot the 95% credible cline region for the selected model
hzar.plot.fzCline(hzar_Hunte$G1_19692754_myosin18ab$analysis$May$model.selected)

####PLOTS####

#"G1_19692754_myosin18ab"####
####Plot all seasons####
#prepare the legend
#G20_13_May=paste0("May (",round(G20_13_May_center,1),", ",round(G20_13_May_width,1),")",", ",round(G20_13_May_AIC))
June=paste0("June (",round(hzar_Hunte$G1_19692754_myosin18ab$analysis$June$modeldetails$param.all$center,1),", ",round(hzar_Hunte$G1_19692754_myosin18ab$analysis$June$modeldetails$param.all$width,1),")")
July=paste0("July (",round(hzar_Hunte$G1_19692754_myosin18ab$analysis$July$modeldetails$param.all$center,1),", ",round(hzar_Hunte$G1_19692754_myosin18ab$analysis$July$modeldetails$param.all$width,1),")")
October=paste0("October (",round(hzar_Hunte$G1_19692754_myosin18ab$analysis$October$modeldetails$param.all$center,1),", ",round(hzar_Hunte$G1_19692754_myosin18ab$analysis$October$modeldetails$param.all$width,1),")")

#plotting
dev.off()
plot.new()
par(mar=c(2, 2, 2, 2))
title(main="Chromosome 1",  adj = 0.95, # Title to the right
      line = -2, outer=TRUE)
mtext("Myosin 18ab - muscle skeletal tissue morphogenesis\n (SNP = G1_19692754)",adj = 0.95,line =-4, outer=TRUE)
#create the area of the left-most plot
par(fig=c(0.01,0.26,0.05,0.9),new=TRUE)
#plot parental Saltwater
plot(x)
plot(hzar_May$AFR[hzar_May$AFR>0.7]~c(1,1,1), xlim=c(0.79,1.51), pch=21,bg="#0e8a81",ylim=c(0,1.05),ylab="", xlab="",yaxt="n", xaxt="n", las=1,add=T)
#Color the panel on the left
rect(0.79,-0.05,1.51,1.09, col = "#70879950") # Color
par(new = TRUE)
plot(hzar_parental_2$AFR[hzar_parental_2$AFR>0.5]~c(1.2,1.2,1.2,1.2,1.2,1.2), xlim=c(0.79,1.51), pch=4, col="grey20", ylim=c(0,1.05),ylab="", xlab="",yaxt="n", xaxt="n", las=1,add=T)

#add axis and tick marks
axis(4,at = c(0, 0.2, 0.4, 0.6, 0.8,1), tck=0.2, las=1, hadj=0.6, cex.axis=0.8)
minor.tick(ny=3, nx=0, tick.ratio=0.5, y.args = list(pos=1.51, lwd=1,col = 1))
mtext(expression(paste( plain("    Frequency\n  Marine allele"))), side=2,line=2.1,padj=2.4,at=0.5,cex=1.2)
par(fig=c(0.159,0.90,0.05,0.9), bg="transparent", new=TRUE)
hzar.plot.fzCline(hzar_Hunte$G1_19692754_myosin18ab$analysis$June$model.selected, type="p",lty=2, pch=15, col="#ffa13a96", fzCol = "#ffffff00", xlim=c(200,300), ylim=c(0,1.05),cex.lab = 1.5, cex.axis = 1, ylab="", xlab="", yaxt="n",xaxt="n")
#hzar.plot.fzCline(G20_13_JunemodelData, type = "p",lty=2, pch = 15,col = "#ffa13a96", fzCol = "#ffffff00", xlim=c(200,300), ylim=c(0,1.05),cex.lab = 1.5, cex.axis = 1, ylab="", xlab="", yaxt="n",xaxt="n")
axis(4, las=1,  tck=0.02, hadj=0.8, cex.axis=0.8)
minor.tick(ny=3, nx=0, tick.ratio=0.5, y.args = list(pos=304, lwd=1,col = 1))
minor.tick(nx=3,ny=3, tick.ratio=0.5)
axis(1, las=1,  tck=0.02, padj=-0.9, cex.axis=0.8)
mtext(expression(paste( plain("Distance (km)"))), side=1,padj=2.2,cex=1.2)

hzar.plot.fzCline(hzar_Hunte$G1_19692754_myosin18ab$analysis$July$model.selected, type="p",lty=2, pch=18, col="#b7410e96", fzCol = "#ffffff00", add=T)
hzar.plot.fzCline(hzar_Hunte$G1_19692754_myosin18ab$analysis$July$model.selected, type="p",lty=2, pch=18, col="#b7410e96", fzCol = "#ffffff00", add=T)

hzar.plot.fzCline(hzar_Hunte$G1_19692754_myosin18ab$analysis$October$model.selected, type="p",lty=2, pch=19, col="#8c8266", fzCol = "#ffffff00", add=T )

legend(x=245,y=1, legend =c(June,July,October), pch=c(15,18,16),lty=c(1,1,1),lwd = c(3,3,3), col = c("#ffa13a", "#b7410e", "#8d8468"), cex = 0.55, title = "Month (center, width)",title.adj = 0.95, pt.cex=1)
#par(fig=c(0,0.8,0,0.8), new=TRUE)
par(fig=c(0.791,1,0.05,0.9),new=TRUE)
plot(hzar_May$AFR[hzar_May$AFR<0.7]~c(1,1,1,1,1), xlim=c(0.99,1.01), pch=21,bg="#855c0b",ylim=c(0,1.05),ylab="", xlab="",yaxt="n", xaxt="n", las=1,add=T)
rect(0.99,-0.05,1.01,1.09, col = "#70879950") # Color
par(new = T)
plot(hzar_parental_3$AFR[hzar_parental_3$AFR>=0]~c(1,1,1,1), xlim=c(0.99,1.01), pch=4, col="grey20", ylim=c(0,1.05),ylab="", xlab="",yaxt="n", xaxt="n", las=1,add=T)
minor.tick(ny=3, nx=0, tick.ratio=0.5, y.args = list(pos=0.99, lwd=1,col = 1))

#"G1_21704763_atp1a1a"####        
#"G2_405619_mucin5"####            
#"G2_10400609_tead1db"####        
#"G2_14538250"####     
#"G4_12806385_EDA"####            
#"G4_19879931"####         
#"G7_1829245"####            
#"G7_20827694_mtnr1bb"####         
#"G8_1302871"####    
#"G9_13197858"####              
#"G10_4238059"####
#"G12_11155054_arhgef25"####   
#"G13_8453167_EPX"####  
#"G14_15149501_SPTLC1"####         
#"G15_7634626"####    
#"G17_6663596_close_to_kif21b"####
#"G18_10525993_zc3h14"####
#"G20_13028350"####











#### corrplot FST####

NA,0.059705286947175,0.12750967243836,0.099145749890343,0.095408844019467,0.161952953564354,0.149770873582322,0.163890514957236,0.127487954307384,0.139960743891531,0.148227053105098,0.059705286947175,NA,0.129560310018232,0.101435282380097,0.088758256097599,0.161589376266388,0.1505102170454,0.163660896490251,0.128340857851225,0.140314145827829,0.148698466219869,0.12750967243836,0.129560310018232,NA,0.094126276337762,0.117970756142902,0.161870137048877,0.144934675247474,0.157313782454421,0.121189245959387,0.136545070083044,0.14271657236478,0.099145749890343,0.101435282380097,0.094126276337762,NA,0.086842099301514,0.143688967765791,0.125707845059751,0.141377919220561,0.095201509132356,0.10973597523623,0.118208985019221,0.095408844019467,0.088758256097599,0.117970756142902,0.086842099301514,NA,0.121516138853939,0.136782265021553,0.151815988067587,0.117351400210947,0.1282700979024,0.13741916346185,0.161952953564354,0.161589376266388,0.161870137048877,0.143688967765791,0.121516138853939,NA,0.171529003722538,0.17971614012259,0.160648223729286,0.163808221300756,0.167928496087163,0.149770873582322,0.1505102170454,0.144934675247474,0.125707845059751,0.136782265021553,0.171529003722538,NA,0.056835827712063,0.14522938739411,0.152445793023959,0.155564596193369,0.163890514957236,0.163660896490251,0.157313782454421,0.141377919220561,0.151815988067587,0.17971614012259,0.056835827712063,NA,0.156771663760888,0.161921616458911,0.165030286856574,0.127487954307384,0.128340857851225,0.121189245959387,0.095201509132356,0.117351400210947,0.160648223729286,0.14522938739411,0.156771663760888,NA,0.099050667483415,0.126637013523647,0.139960743891531,0.140314145827829,0.136545070083044,0.10973597523623,0.1282700979024,0.163808221300756,0.152445793023959,0.161921616458911,0.099050667483415,NA,0.123962403085312,0.148227053105098,0.148698466219869,0.14271657236478,0.118208985019221,0.13741916346185,0.167928496087163,0.155564596193369,0.165030286856574,0.126637013523647,0.123962403085312,NA

FST_pwc<-structure(c(NA,0.059705286947175,0.12750967243836,0.099145749890343,0.095408844019467,0.161952953564354,0.149770873582322,0.163890514957236,0.127487954307384,0.139960743891531,0.148227053105098,0.059705286947175,NA,0.129560310018232,0.101435282380097,0.088758256097599,0.161589376266388,0.1505102170454,0.163660896490251,0.128340857851225,0.140314145827829,0.148698466219869,0.12750967243836,0.129560310018232,NA,0.094126276337762,0.117970756142902,0.161870137048877,0.144934675247474,0.157313782454421,0.121189245959387,0.136545070083044,0.14271657236478,0.099145749890343,0.101435282380097,0.094126276337762,NA,0.086842099301514,0.143688967765791,0.125707845059751,0.141377919220561,0.095201509132356,0.10973597523623,0.118208985019221,0.095408844019467,0.088758256097599,0.117970756142902,0.086842099301514,NA,0.121516138853939,0.136782265021553,0.151815988067587,0.117351400210947,0.1282700979024,0.13741916346185,0.161952953564354,0.161589376266388,0.161870137048877,0.143688967765791,0.121516138853939,NA,0.171529003722538,0.17971614012259,0.160648223729286,0.163808221300756,0.167928496087163,0.149770873582322,0.1505102170454,0.144934675247474,0.125707845059751,0.136782265021553,0.171529003722538,NA,0.056835827712063,0.14522938739411,0.152445793023959,0.155564596193369,0.163890514957236,0.163660896490251,0.157313782454421,0.141377919220561,0.151815988067587,0.17971614012259,0.056835827712063,NA,0.156771663760888,0.161921616458911,0.165030286856574,0.127487954307384,0.128340857851225,0.121189245959387,0.095201509132356,0.117351400210947,0.160648223729286,0.14522938739411,0.156771663760888,NA,0.099050667483415,0.126637013523647,0.139960743891531,0.140314145827829,0.136545070083044,0.10973597523623,0.1282700979024,0.163808221300756,0.152445793023959,0.161921616458911,0.099050667483415,NA,0.123962403085312,0.148227053105098,0.148698466219869,0.14271657236478,0.118208985019221,0.13741916346185,0.167928496087163,0.155564596193369,0.165030286856574,0.126637013523647,0.123962403085312,NA), .Dim=c(11,11), .Dimnames=list(c("Ems","Weser","HGLA","LM","HB","RM","HBU","HBO","BrB","AM","FB"),c("Ems","Weser","HGLA","LM","HB","RM","HBU","HBO","BrB","AM","FB")))

FST_pwc<-structure(c(NA,0.1497708736,0.05970528695,0.09540884402,0.1275096724,0.1399607439,0.1274879543,0.1482270531,0.163890515,0.1619529536,0.09914574989,0.1497708736,NA,0.150510217,0.136782265,0.1449346752,0.152445793,0.1452293874,0.1555645962,0.05683582771,0.1715290037,0.1257078451,0.05970528695,0.150510217,NA,0.0887582561,0.12956031,0.1403141458,0.1283408579,0.1486984662,0.1636608965,0.1615893763,0.1014352824,0.09540884402,0.136782265,0.0887582561,NA,0.1179707561,0.1282700979,0.1173514002,0.1374191635,0.1518159881,0.1215161389,0.0868420993,0.1275096724,0.1449346752,0.12956031,0.1179707561,NA,0.1365450701,0.121189246,0.1427165724,0.1573137825,0.161870137,0.09412627634,0.1399607439,0.152445793,0.1403141458,0.1282700979,0.1365450701,NA,0.09905066748,0.1239624031,0.1619216165,0.1638082213,0.1097359752,0.1274879543,0.1452293874,0.1283408579,0.1173514002,0.121189246,0.09905066748,NA,0.1266370135,0.1567716638,0.1606482237,0.09520150913,0.1482270531,0.1555645962,0.1486984662,0.1374191635,0.1427165724,0.1239624031,0.1266370135,NA,0.1650302869,0.1679284961,0.118208985,0.163890515,0.05683582771,0.1636608965,0.1518159881,0.1573137825,0.1619216165,0.1567716638,0.1650302869,NA,0.1797161401,0.1413779192,0.1619529536,0.1715290037,0.1615893763,0.1215161389,0.161870137,0.1638082213,0.1606482237,0.1679284961,0.1797161401,NA,0.1436889678,0.09914574989,0.1257078451,0.1014352824,0.0868420993,0.09412627634,0.1097359752,0.09520150913,0.118208985,0.1413779192,0.1436889678,0), .Dim=c(11,11), .Dimnames=list(c("Ems","HBU","Weser","HB","LM","AM","BrB","FB","HBO","RM","HGLA"),c("Ems","HBU","Weser","HB","LM","AM","BrB","FB","HBO","RM","HGLA")))
low <- matrix(NA, nrow = nrow(FST_pwc), ncol = ncol(FST_pwc))
up <- matrix(NA, nrow = nrow(FST_pwc), ncol = ncol(FST_pwc))

#populate upper and lower matrices
up[upper.tri(up)] <- FST_pwc[upper.tri(FST_pwc)]
low[lower.tri(low)] <- FST_pwc[lower.tri(FST_pwc)]

#pivot upper and lower for plotting
lower_dat <- low|>
  as.data.frame() |>
  `colnames<-`(colnames(FST_pwc)) |>
  mutate(xvar = colnames(FST_pwc)) |>
  pivot_longer(cols = -xvar, names_to = "yvar") 

upper_dat <- up|>
  as.data.frame() |>
  `colnames<-`(colnames(FST_pwc)) |>
  mutate(xvar = colnames(FST_pwc)) |>
  pivot_longer(cols = -xvar, names_to = "yvar") 

lower_dat|> #lower matrix data
  ggplot(aes((xvar), yvar))+ 
  #geom_tile(fill = NA, color = "grey")+ #background grid
  #geom_point(aes(fill = value, size = value), pch = 22)+ # differnt sized points
  geom_text(data = upper_dat, aes(color = value, label = round(value, 2)))+ #plot cor in upper right
  scale_size_continuous(breaks = seq(-1, 1, by = 0.5))+ # define size breaks
  labs(x = "", y = "")+ #remove unnecessary labels
  #scale_fill_gradient2(low = "darkred",mid = "white", high = "darkblue", midpoint = 0)+ #define square colors
  scale_color_gradient2(low = "darkred",mid = "white", high = "darkblue", midpoint = 0)+ #define text colors
  scale_x_discrete(limits = rev)+# rev to make the triagle a certain side
  #make it look pretty
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))

dist2order = function(corr, method, ...) {
  d_corr = as.dist(1 - corr)
  s = seriate(d_corr, method = method, ...)
  i = get_order(s)
  return(i)
}

corrplot(FST_pwc, type="upper", is.corr = FALSE, order=,col = colorRampPalette(c("grey75","grey85","grey95","#d8f3dc","#b7e4c7","#95d5b2","#74c69d","#52b788","#40916c","#2d6a4f"))(100),
         addCoef.col = TRUE,method = "shade",tl.col = "black",diag=F, title="Genome-wide mean FST values for PoolSeq",mar=c(0,0,1,0) ) %>% corrRect(name=c("HB","RM","HBU","HBO","BrB","AM","FB"))
i=dist2order (FST_pwc, method="Random")
corrplot(FST_pwc[i,i],cl.pos=n)

i=c("Ems","Weser","HGLA","LM","HB","RM","HBU","HBO","BrB","AM","FB")
c("grey80","#ffba08","#faa307","#f48c06","#e85d04","#dc2f02","#d00000","#9d0208","#6a040f","#370617")
corrplot(FST_pwc, method="square", type = "lower", tl.pos = "l", cl.pos = "n", add = TRUE)

0,0.1497708736,0.05970528695,0.09540884402,0.1275096724,0.1399607439,0.1274879543,0.1482270531,0.163890515,0.1619529536,0.09914574989,0.1497708736,0,0.150510217,0.136782265,0.1449346752,0.152445793,0.1452293874,0.1555645962,0.05683582771,0.1715290037,0.1257078451,0.05970528695,0.150510217,0,0.0887582561,0.12956031,0.1403141458,0.1283408579,0.1486984662,0.1636608965,0.1615893763,0.1014352824,0.09540884402,0.136782265,0.0887582561,0,0.1179707561,0.1282700979,0.1173514002,0.1374191635,0.1518159881,0.1215161389,0.0868420993,0.1275096724,0.1449346752,0.12956031,0.1179707561,0,0.1365450701,0.121189246,0.1427165724,0.1573137825,0.161870137,0.09412627634,0.1399607439,0.152445793,0.1403141458,0.1282700979,0.1365450701,0,0.09905066748,0.1239624031,0.1619216165,0.1638082213,0.1097359752,0.1274879543,0.1452293874,0.1283408579,0.1173514002,0.121189246,0.09905066748,0,0.1266370135,0.1567716638,0.1606482237,0.09520150913,0.1482270531,0.1555645962,0.1486984662,0.1374191635,0.1427165724,0.1239624031,0.1266370135,0,0.1650302869,0.1679284961,0.118208985,0.163890515,0.05683582771,0.1636608965,0.1518159881,0.1573137825,0.1619216165,0.1567716638,0.1650302869,0,0.1797161401,0.1413779192,0.1619529536,0.1715290037,0.1615893763,0.1215161389,0.161870137,0.1638082213,0.1606482237,0.1679284961,0.1797161401,0,0.1436889678,0.09914574989,0.1257078451,0.1014352824,0.0868420993,0.09412627634,0.1097359752,0.09520150913,0.118208985,0.1413779192,0.1436889678,0
                   
