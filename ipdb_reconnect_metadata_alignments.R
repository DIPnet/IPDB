############################
# A script for taking metadata from cleaned and QC'd metadata files (mdfastaQC2) 
# and joining it with cleaned, QC'd and aligned mtDNA data (alignedQC2). 
#The script jettisons the original, non-aligned mtDNA data that is in the mdfastaQC2 files.
# By Eric Crandall, Cynthia Riginos & Libby Liggins. January 2014 through May 2015.
############################


library(plyr)
datapath<-"~/Google Drive/!DIPnet_DB/Repository"
localpath<-"~/Desktop/Reunite_metadata_and_alignments" #or whatever
setwd(path)

#new folder for downstream db
downstream<-file.path(localpath,Sys.Date())
dir.create(downstream)

#here are the columns we'll be keeping if we need a smaller file
keeps<-c("IPDB_ID", "Genus_species_locus","materialSampleID","genus","species","phylum","locus","decimalLatitude","decimalLongitude","georeferenceProtocol","country","island","stateProvince","locality","recordedBy","principalInvestigator","yearCollected","sequence") 

##LOAD IN QCd METADATA ##
#initialize an empty data frame with as many columns as are in keeps
allmetadataframe<-data.frame(matrix(ncol=0, nrow=0))

filecount<-0

#loop through all the files and bung'em together
for(file in (list.files(paste(datapath,"/1-cleaned_QC2_metadata_files", sep=""),pattern="mdfasta",full.names=F))){ #file<-"name_checked_mdfasta_trimax_CO1_td.txt"
    
  print(file)
  data<-read.table(file=paste("1-cleaned_QC2_metadata_files", file, sep="/"),header=T,stringsAsFactors=F,sep="\t",na.strings=c("","NA","#N/A"),strip.white=T,fill=T,comment.char="",quote="", colClasses=c("materialSampleID"="character"))

  #break up the file name and pull out the locus name, convert it to upper case
	locus<-strsplit(file,split="_")[[1]][3]
	locus<-toupper(locus)
  locus<-gsub("COI", "CO1", locus)  
  locus<-gsub("ATPASE", "A68", locus)  
  locus<-gsub("CYTB", "CYB", locus)  
  
	#add the locus (and file name while we are troubleshooting)
	data<-cbind(data[,1:(length(data)-1)],locus,file,sequence=data[,length(data)])
		
  #bung it on there!

  allmetadataframe<-join(allmetadataframe, data, type = "full")
  filecount<-filecount+1
		
}
print(paste(filecount,"files processed"))


##LOAD IN QC'd ALIGNMENT DATA ##
allalignframe<-data.frame(matrix(ncol=0, nrow=0))

for(fasfile in (list.files(paste(datapath,"/1-cleaned_fasta_files", sep=""),pattern="aligned",full.names=F))){ 
  print(fasfile)
  fas<-readLines(paste("1-cleaned_fasta_files", fasfile, sep="/"))

  #break up the file name and pull out the species, genus, and locus names
  genus<-strsplit(fasfile,split="_")[[1]][2]
  species<-strsplit(fasfile,split="_")[[1]][3]
  locus<-strsplit(fasfile,split="_")[[1]][4]
  locus<-gsub("(\\w+)\\..*", "\\1", locus)  #removes ".fas" or ".fasta" etc
  locus<-toupper(locus)
  locus<-gsub("COI", "CO1", locus)  
  locus<-gsub("ATPASE", "A68", locus)  
  locus<-gsub("CYTB", "CYB", locus)  
  
  
  #grab fasta IDS as odd lines and sequences as even lines make fasframe dataframe
  names <- fas[seq(1, length(fas), 2)]
  seqs<- fas[seq(2, length(fas), 2)]
  materialSampleID<-gsub(paste(">(.*)_", genus, ".*", sep=""), "\\1", names)  #pulls out front of name when genus is included
  materialSampleID<-gsub(">(.*)", "\\1", materialSampleID)   #if genus is not included, then the ">" remains and needs to be removed
  fasframe<-data.frame(materialSampleID)
  fasframe$materialSampleID<-as.character(materialSampleID)
  fasframe$nameInAlignmentFile<-gsub(">(.*)", "\\1", names) 
  fasframe$alignmentfile<-fasfile
  fasframe$IPDB_ID<-paste(genus, species, locus,fasframe[,"materialSampleID"], sep="_")
  fasframe$Genus_species_locus<-paste(genus, species, locus, sep="_")
  
  fasframe$sequence<-seqs
  
  #join to alignment frame and tidy up
  allalignframe<-join(allalignframe, fasframe, type = "full")
    
}

##OPTIONAL - CREATE FRAMES TO ADD IN SPECIAL ALIGNMENTS###
#specialtotalframe<-data.frame(matrix(ncol=0, nrow=0))
#for(file in (list.files(paste(path,"/Special_species/completed", sep=""),pattern="totaldata",full.names=F))){ 
  
#  print(file)
#  specialdata<-read.table(file=paste("Special_cases_completed_reuniting", file, sep="/"),header=T,stringsAsFactors=F,sep="\t",na.strings=c("","NA","#N/A"),strip.white=T,fill=T,comment.char="",quote="")
#  specialdata[,"Genus_species_locus"]<-NULL
#  specialtotalframe<-join(specialtotalframe, specialdata, type="full")
#}
#specialtotalframe$Genus_species_locus<-paste(specialtotalframe[,"genus"],specialtotalframe[,"species"],specialtotalframe[,"locus"], sep="_")
#write.table(specialtotalframe,file=paste(path,"/",Sys.Date(),"/",Sys.Date(),"_specialtotalframe.txt",sep=""),quote=F,sep="\t",row.names=F)


specialmatchedframe<-data.frame(matrix(ncol=0, nrow=0))
for(file in (list.files(paste(path,"/Special_species/completed", sep=""),pattern="matcheddata",full.names=F))){ 
  
  print(file)
  specialdata<-read.table(file=paste(path,"Special_species/completed",file, sep="/"),header=T,stringsAsFactors=F,sep="\t",na.strings=c("","NA","#N/A"),strip.white=T,fill=T,comment.char="",quote="")
  specialdata[,"Genus_species_locus"]<-NULL
  specialmatchedframe<-join(specialmatchedframe, specialdata, type = "full")
}


specialmatchedframe$Genus_species_locus<-paste(specialmatchedframe[,"genus"],specialmatchedframe[,"species"],specialmatchedframe[,"locus"], sep="_")
#write.table(specialmatchedframe,file=paste(path,"/",Sys.Date(),"/",Sys.Date(),"_specialmatchedframe.txt",sep=""),quote=F,sep="\t",row.names=F)



##COMBINE METADATA AND ALIGNMENTS##

#add a unique identifier based on matID+genus+species+locus to metadata dataframe
#allmetadataframe$materialSampleID<-as.character(allmetadataframe$materialSampleID) #this is now taken care of in the read.table() statement
IPDB_ID<-paste(allmetadataframe[,"genus"],allmetadataframe[,"species"],allmetadataframe[,"locus"],allmetadataframe[,"materialSampleID"], sep="_")
Genus_species_locus<-paste(allmetadataframe[,"genus"],allmetadataframe[,"species"],allmetadataframe[,"locus"], sep="_")
allmetadataframe<-cbind(IPDB_ID,Genus_species_locus,allmetadataframe)
allmetadataframe[,"sequence"]<-NULL   #sequences removed because they will be replaced with aligned sequences

#combine files 
totaldataframe<-join(allmetadataframe, allalignframe, type = "full")
matcheddataframe<-join(allmetadataframe, allalignframe, type = "inner")

#if adding special alignments
totaldataframe<-join(totaldataframe, specialmatchedframe, type = "full")
matcheddataframe<-join(matcheddataframe, specialmatchedframe, type = "full")


keepsdataframe<-matcheddataframe[,keeps]

#write to a time stamped directory
write.table(allmetadataframe,file=paste(localpath,"/",Sys.Date(),"/",Sys.Date(),"_allmetadata.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(allalignframe,file=paste(localpath,"/",Sys.Date(),"/",Sys.Date(),"_allalignmentdata.txt",sep=""),quote=F,sep="\t",row.names=F)

write.table(totaldataframe,file=paste(localpath,"/",Sys.Date(),"/",Sys.Date(),"_totaldata.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(matcheddataframe,file=paste(localpath,"/",Sys.Date(),"/",Sys.Date(),"_matcheddata.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(keepsdataframe,file=paste(localpath,"/",Sys.Date(),"/",Sys.Date(),"_keepsdata.txt",sep=""),quote=F,sep="\t",row.names=F)

##IDENTIFYING TROUBLE FILES##
##species with TRUE values indicate that either metadata or alignments are missing
#missing sequences
with(totaldataframe, table(Genus_species_locus,(is.na(sequence))))
#missing metadata
with(totaldataframe, table(Genus_species_locus,(is.na(principalInvestigator))))


#find sequences missing from a particular species
totaldataframe[totaldataframe$Genus_species_locus=="Labroides_dimidiatus_CR" & is.na(totaldataframe$sequence),	3]

#find metadata missing from a particular species
totaldataframe[totaldataframe$Genus_species_locus=="Zebrasoma_flavescens_CYB" & is.na(totaldataframe$principalInvestigator),  7]

# The above script requires that fasta files contain samples from only a single species, so here is a script to break multi-species fasta files into their component species

#load in each multispecies pair of files (not going to loop this)
multimeta1<-read.table(file="/Users/eric/Google Drive/!DIPnet_DataQC/Reunite_metadata_and_alignments/1-cleaned_QC2_metadata_files/mdfastaQC2_ptespp_16S_MK.txt",sep="\t",header=T,stringsAsFactors = F)
multifas1<-read.fasta(file="/Users/eric/Google Drive/!DIPnet_DataQC/Reunite_metadata_and_alignments/Put missing files HERE/multispecies/multispecies_ptespp_16S_MK.fasta",forceDNAtolower = F)

#find all unique species
species<-unique(paste(multimeta1$genus,multimeta1$species))

#pick out the msids for each species from the metadata then use them to subset the fasta object and write to a file
for(s in species){
  msids<-multimeta1$materialSampleID[paste(multimeta1$genus,multimeta1$species)==s]
  subfas<-multifas1[msids]
  write.fasta(subfas,names=names(subfas),file.out=paste(paste("alignedQC",strsplit(s," ")[[1]][1],strsplit(s," ")[[1]][2],"16S",sep="_"),"fasta",sep="."),nbchar=1000)
}
