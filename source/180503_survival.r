cd("/Users/briglo/Desktop/briglo_backup/scRNA/davgal/paperstuff/r_objects")
library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")
# Test the CGDS endpoint URL using a few simple API tests
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)[,1]


# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[30,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[4,1]

# Get data slices for a specified list of genes, genetic profile and case list

load("180502_epires1.2posmarkers.rdata")
dat<-lapply(smark,function(x) getProfileData(mycgds,x,mygeneticprofile,mycaselist))
clin<-getClinicalData(mycgds,mycaselist)
save(smark,dat,clin,file="180502_cgdsrData.rdata")


pcas<-lapply(dat,function(x) {
tx<-x[rowSums(!is.na(x))>0,colSums(!is.na(x))>0]
tx[is.na(tx)]=0
return(prcomp(tx))
}
)

par(mfrow=c(3,5))
lapply(pcas,function(x) plot(density(x$x[,1])))
dev.copy2pdf(file="180502_epires1.2marker_TCGA_pc1_dist.pdf")

cutoffs<-c(-3,-4,5,-1,-4,-2.5,-1,-4,0,-3,2,0,0,0,-4)
par(mfrow=c(3,5))
lapply(pcas,function(x) plot(x$x[,1],x$x[,2])))


classes<-lapply(1:15,function(i){
	return(data.frame(pcas[[i]]$x[,1]>cutoffs[i]))
})

par(mfrow=c(3,5))
lapply(classes, function(x){
tmp<-clin[rownames(x),]
 plot(survfit(Surv(as.numeric(tmp$OS_MONTHS),ifelse(tmp$OS_STATUS=="DECEASED",1,0))~x[,1]))


 text(0,0,summary(coxph(Surv(as.numeric(tmp$OS_MONTHS),ifelse(tmp$OS_STATUS=="DECEASED",1,0))~x[,1]))$logtest[3],pos=4)
}
)
dev.copy2pdf(file="180502_epires1.2marker_TCGA_survPlitByPCA.pdf")


#trying by all
allmark<-unique(unlist(smark))
tmp<-lapply((0:9*200)+1, function(x) getProfileData(mycgds,allmark[seq(x,x+199,1)],mygeneticprofile,mycaselist))
x<-data.frame(do.call(cbind,tmp))
tx<-x[rowSums(!is.na(x))>0,colSums(!is.na(x))>0]
tx[is.na(tx)]=0
pca_all<-prcomp(tx)
head(pca_all$x)


plot(pca_all$x[,1],pca_all$x[,2]) ; abline(a=-15,b=-.4)
colvec<-ifelse(pca_all$x[,2]>-.3*pca_all$x[,1]-18,"blue",'red')
plot(pca_all$x[,1],pca_all$x[,2],col=colvec) ; abline(a=-18,b=-.3)
tmp<-clin[names(colvec),]
plot(survfit(Surv(as.numeric(tmp$OS_MONTHS),ifelse(tmp$OS_STATUS=="DECEASED",1,0))~colvec))
text(0,0,summary(coxph(Surv(as.numeric(tmp$OS_MONTHS),ifelse(tmp$OS_STATUS=="DECEASED",1,0))~colvec))$logtest[3],pos=4)



