HANDOVER chrcha

1.2 source("source/endTOend.r")

source("source/seuratFunctions.r")
source("source/scGPSfunctions.r")
load("r_objects/cell.cyclegenes.rdata")
source/190613_runscGPS.r

load("190504_magicALLgenes.rdata")
load("r_objects/190524_magicPHATEall.rdata")

load("r_objects/190614_scGPS_byLocation.rdata") #check this is output from runscgps


save(comb,file="BYLOCATION.rdata")
save(prim,mets,file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190614_scGPSinputData.rdata")
save(DEgenes,file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190613_DEgenesByPTcluster")
save(DEgenes,file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190613_DEgenesByLungCluster.rdata")
save(md,file="r_objects/190606_4by4_IntegratedSeurat_metadata.rdata")
save(integrated,file="r_objects/190606_4by4_IntegratedSeurat.rdata")


#To Do chrcha
add magic and scent methods
add AR analysis methods
plot scgps methods
get enrichment function sorted for chrcha object
list of the top genes from the scGPS model for the cluster pairs that have a Pr > 0.5