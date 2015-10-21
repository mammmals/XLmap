#Contact map score (CMScore)
#Script to find matched presence in two matrices/arrays
#Initial purpose for comparing cross-linked protein data to contact maps
#CMS=contact map score
# '#print' are used for troubleshooting various steps of the script

cms.example=function(distancerange,exampleinput) {
	if(is.numeric(distancerange)==FALSE){
		return("Error: Contact map distances must be numeric.")
	}else if(is.character(exampleinput)==FALSE){
		return("Error: Protein input must be character, add quote marks around protein name/accession.")
	}else{
		#libraries
			#library(Matrix)
			#library(zoo)
			#library(Rpdb)
			#library(lattice)
			#library(latticeExtra)
			#library(tcltk)
		#Nullify
		pdbfile.initial=NULL
		pdbfile=NULL
		xatoms=NULL
		xmolecule=NULL
		xchaina=NULL
		xalpha=NULL
		aamin=NULL
		xlstart=NULL
		xlend=NULL
		xdistmatrix=NULL
		distcon=NULL
		colfunc1=NULL
		colfunc2=NULL
		hmcols=NULL
		xlcontactfile=NULL
		uniprotsitetable=NULL
		upaccs1=NULL
		upaccs2=NULL
		upinputselect1=NULL
		upinputselect2=NULL
		reciprocalust=NULL
		crosslinkmap=NULL
		upaccs1=NULL
		upaccs2=NULL
		xlcontactfile=NULL
		uniprotsitetable=NULL
		upinputselect1=NULL
		upinputselect2=NULL
		reciprocalust=NULL
		crosslinktable=NULL
		crosslinkmap=NULL
		sort.clt=NULL
		cms.crosslinkmap=NULL
		cms.matrix.input=NULL
		cms.matrix=NULL
		cmscore.table=NULL
		
		recname=chainid=elename=NULL
		
		#data
		#data(cmdata)
		#data(xlmapdata)
		
		#import pdb from example file
		pdbfile.initial=paste("p.",exampleinput,sep="")
		pdbfile=get(pdbfile.initial)

		
		cmscore.table.initial=unlist(
			lapply(distancerange,function(cadistance){
			
			###
			###Making PDB matrix
			##
			##if pdb file imported, find contacts within x angstroms and plot them to data matrix
			#first identify atoms
			xatoms=pdbfile$atoms
			
			#as.data.frame important for maintaining values in eleid and resid columns
			xmolecule=subset(as.data.frame(xatoms),recname=="ATOM")
			
			yes.pdb=c('oxa23','tbpa','ppib')
			no.pdb=c('ab57.2983','ppib.decoy')
			
			#then identify alpha-carbons in chain A
			if(any(exampleinput==yes.pdb)==TRUE) {
				xchaina=subset(xmolecule,chainid=="A")
				xalpha=subset(xchaina,elename=="CA")
			}else if(any(exampleinput==no.pdb)==TRUE){
				xalpha=subset(xmolecule,elename=="CA")
			}else {
			return("No example files with that name.")
			}
						
				
			aamin=1
			xlstart=min(xalpha[c("resid")])
			xlend=max(xalpha[c("resid")])
			
			#create matrix including any crystal structure offset in amino acid sequence
			xfill=as.data.frame(c(1:xlend))
			colnames(xfill)=c('resid')
			xalphafilled.initial=merge(xfill,xalpha,all=TRUE)
			xalphafilled=xalphafilled.initial[1:xlend,]
			
			#determine distance between alpha carbons (this will create the matrix of interactions)
			xdistmatrix=as.matrix(dist(xalphafilled[c("x1","x2","x3")]))
			xdistmax=max(as.matrix(dist(xalpha[c("x1","x2","x3")])))
			if(cadistance>max(xdistmax)){
				cadistance=max(xdistmax)
			}else{
				cadistance=cadistance
			}
			xdistmatadj=xdistmatrix[xdistmatrix>cadistance]=50
			
			cms.xdm=xdistmatrix
			cms.pdb.distmatadj=cms.xdm[cms.xdm<cadistance]=1

			
				
			###
			###
			###
			###Make matrix for XL data from example data
			###
			###
			
			crosslinktable=get(gsub('.decoy',"",exampleinput))
			crosslinkmap=as.matrix(crosslinktable[c(3,4)])

			
			#Make Xl matrix for ContactMapScore(CMS)
			sort.clt=crosslinktable[order(crosslinktable$mod_pos1,crosslinktable$mod_pos2),]
			cms.crosslinkmap=as.matrix(sort.clt[c(3,4)])
			cms.matrix.initial=as.data.frame(cms.crosslinkmap)
			#remove crosslinks greater than xlend
			cms.matrix.namaker=cms.matrix.initial[cms.matrix.initial>xlend]=NA
			cms.matrix.short=na.omit(cms.matrix.initial)
			#fill to max residue in pdb file
			cms.matrix.input=rbind(cms.matrix.short,c(xlend,xlend))

			cms.matrix=Matrix(0,nrow=nrow(cms.xdm),ncol=ncol(cms.xdm))
			cms.matrix[cbind(cms.matrix.short$mod_pos1,cms.matrix.short$mod_pos2)]=1
			cms.matrix=as.matrix(cms.matrix)
			
			#
			###
			###
			###
			#Comparing the two matrices
			###
			###
			###
			#
				
			
			#change row and column names in both matrices to match each other
			rownames(cms.xdm)=colnames(cms.xdm)=rownames(cms.matrix)=colnames(cms.matrix)=c(1:xlend)
			
			#find matrix sum
			matrix.matcher=cms.xdm==cms.matrix
			matrix.matcher=matrix.matcher+0
			matrix.overlap=sum(matrix.matcher,na.rm=TRUE)
			
			cms.crosslink.sum=sum(cms.matrix)
			
			cms.xdm.sum.input=cms.xdm
			cms.xdm.sum.middle=cms.xdm.sum.input[cms.xdm.sum.input!=1]=0
			cms.xdm.sum=sum(cms.xdm.sum.input,na.rm=TRUE)
			
			cms.middle=cms.crosslink.sum/matrix.overlap
			#Lower scores are better, 0 is perfect overlap!
			cms.final=cms.middle-1
			
			return(cms.final)
		}
		)
		)
		#create plot and add CMScore as text to plot, users can save this as a metafile to get necessary information;
		if(is.infinite(max(cmscore.table.initial))) {
			cmscore.table.maxfind=cmscore.table.initial
			cmscore.table.maxfind[is.infinite(cmscore.table.maxfind)]=0
			max.cmscore.adj=max(cmscore.table.maxfind)
			cmscore.table.final=cmscore.table.initial
			cmscore.table.final[is.infinite(cmscore.table.final)]=max.cmscore.adj
		} else {
			cmscore.table.final=cmscore.table.initial
		}
	#Plot the residual scores and CMScore
		cmscore.table.final=cbind(c(min(distancerange):max(distancerange)),cmscore.table.final)
		plot(
			cmscore.table.final,type='o',
			xlab="Distance Constraint used to Generate Contact Map (Angstroms)",
			ylab="Residual Score",
			xlim=c(min(distancerange),max(distancerange)),
			main=paste("CMScore =",round(sum(cmscore.table.final[,2]),digits=5))
		)
		question.save.plot = winDialog(type = c("yesno"),"Do you want to save the CMScore table?")
		if(question.save.plot=="YES") {
			write.table(cmscore.table.final,file = paste(exampleinput,"CMScoreTable",format(Sys.time(),"%Y%m%d%H%M%S"),".txt", sep=""), sep="\t",row.names=F)
			print("NICE PLOT!")
		} else {
			print("NICE PLOT!")
		}
	}
}

#cmsexample(0:40,'ppib')
