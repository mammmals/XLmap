cmscore=function(distancerange,cms.proteininput, molfile=NULL, xltable=NULL,
					cms.pdbx=FALSE,cms.offset=0,cms.sites.only=FALSE,cms.which.chain=NULL, 
					save.cms.table=TRUE, ...){
	if(is.numeric(distancerange)==FALSE){
		return("Error: Contact map distances must be numeric.")
	}else if(is.character(cms.proteininput)==FALSE){
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
		
		if(is.null(cms.which.chain)) {
			cms.which.chain="A"
		} else {
			cms.which.chain=cms.which.chain
		}
		
		
		#import pdb or cms.pdbx file from user chosen directory
		if(cms.pdbx==FALSE) {
			if(is.null(molfile)){
				pdbfile=read.pdb(tk_choose.files(caption="Choose a .pdb file to create a contact map."))
			} else {
				pdbfile=read.pdb(molfile)
			}
			#first identify atoms
			xatoms=pdbfile$atoms
			#as.data.frame important for maintaining values in eleid and resid columns
			xmolecule=subset(as.data.frame(xatoms),recname=="ATOM")
		} else {
			if(is.null(molfile)){
				pdbx.atoms=read.pdbx.atoms(tk_choose.files(caption="Choose PDBx/mmCIF file for contact map."))
			} else {
				pdbfile=read.pdbx.atoms(molfile)
			}
			colnames(cms.pdbx.atoms)=c(
				'recname',#'group_PDB',
				'id',
				'type_symbol',
				'elename',#'label_atom_id',
				'label_alt_id',
				'label_comp_id',
				'chainid',#'label_asym_id',
				'label_entity_id',
				'label_seq_id',
				'cms.pdbx_PDB_ins_code',
				'x1',#'Cartn_x',
				'x2',#'Cartn_y',
				'x3',#'Cartn_z',
				'occupancy',
				'B_iso_or_equiv',
				'Cartn_x_esd',
				'Cartn_y_esd',
				'Cartn_z_esd',
				'occupancy_esd',
				'B_iso_or_equiv_esd',
				'cms.pdbx_formal_charge',
				'resid',#'auth_seq_id',
				'auth_comp_id',
				'auth_asym_id',
				'auth_atom_id',
				'cms.pdbx_PDB_model_num'
			)
			xmolecule=cms.pdbx.atoms
		}
		
		#then identify alpha-carbons in the user defined chain
		xchaina=subset(xmolecule,chainid==cms.which.chain)
		xalpha=subset(xchaina,elename=="CA")
		
		aamin=1
		xalpha.resid=xalpha[c("resid")]
		xalpha.resid[]=lapply(xalpha.resid,as.numeric)
		xlstart=min(xalpha.resid)
		xlend=max(xalpha.resid)
		
		#create matrix including any crystal structure offset in amino acid sequence
		xfill=as.data.frame(c(1:xlend))
		colnames(xfill)=c('resid')
		xalphafilled.initial=merge(xfill,xalpha,all=TRUE)
		xalphafilled=xalphafilled.initial[1:xlend,]
			
		#input is uniprot ID
		if(is.null(xltable)){
			xlcontactfile=read.delim(
				tk_choose.files(caption="Choose text file containing crosslinked relationships."),as.is=TRUE)
		} else {
			xlcontactfile=read.delim(xltable,as.is=TRUE)
		}
		
		#count total number of unique crosslinks
		count.links=xlcontactfile[order(xlcontactfile$mod_pos1,xlcontactfile$mod_pos2),]
		count.links=unique(count.links)
		total.crosslinks=nrow(count.links)
		
		cmscore.table.initial=unlist(
			lapply(distancerange,function(cadistance){
			
			###
			###Making PDB matrix
			##
			#determine distance between alpha carbons (this will create the matrix of interactions)
			xdistmatrix=as.matrix(dist(xalphafilled[c("x1","x2","x3")]))
			xdistmax=max(as.matrix(dist(xalpha[c("x1","x2","x3")])))
			if(cadistance>max(xdistmax)){
				cadistance=max(xdistmax)
				#print(cadistance)
				#print(max(xdistmax))
			}else{
				cadistance=cadistance
				#print(cadistance)
				#print(max(xdistmax))
			}
			xdistmatadj=xdistmatrix[xdistmatrix>cadistance]=50
			
			cms.xdm=xdistmatrix
			cms.pdb.distmatadj=cms.xdm[cms.xdm<cadistance]=1

			###
			###
			###
			###Make matrix for XL data from react2 files
			###
			###

			#input as uniprot accession or as sequest output
			#if a prot1 column contains a pipe then run grepregx otherwise take the column as is
			
			if(grepl("|",xlcontactfile[1,c('prot1')],fixed=TRUE)){
				#Find  Uniprot Accession number between pipes and and site of modification then output the results
				upaccs1=sapply(xlcontactfile[c('prot1')],function(y) {
							substring(
								y,
								first = 4,
								last = gregexpr("\\Q|\\E.",y)[[1]][c(2)]-1
							)
					}
				)

				#Find  Uniprot Accession number between pipes and and site of modification then output the results
				upaccs2=sapply(xlcontactfile[c('prot2')],function(y) {
						substring(
							y,
							first = 4,
							last = gregexpr("\\Q|\\E.",y)[[1]][c(2)]-1
						)
					}
				)
			} else {
				upaccs1=xlcontactfile[c('prot1')]
				upaccs2=xlcontactfile[c('prot2')]
			}
			upaccs1=as.vector(upaccs1)
			upaccs2=as.vector(upaccs2)
			

			#Creates intraprotein link table based on cms.proteininput
			#
			# AND based on cms.offset
			#
			xlcontactfile[c('mod_pos1')]=xlcontactfile[c('mod_pos1')]+cms.offset
			xlcontactfile[c('mod_pos2')]=xlcontactfile[c('mod_pos2')]+cms.offset
			
			if(cms.sites.only==FALSE) {
				uniprotsitetable=as.data.frame(
					cbind(upaccs1,upaccs2,xlcontactfile[c('mod_pos1')],xlcontactfile[c('mod_pos2')])
				)
				#Double the subsetting
				#print(cms.proteininput)
				upinputselect1=subset(uniprotsitetable,upaccs1==cms.proteininput)
				upinputselect2=subset(upinputselect1,upaccs2==cms.proteininput)
				reciprocalust=upinputselect2[c(1,2,4,3)]
				colnames(reciprocalust)=colnames(upinputselect2)
				crosslinktable=rbind(upinputselect2,reciprocalust)
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
			} else {
				uniprotsitetable=as.data.frame(cbind(xlcontactfile[c('mod_pos1')],xlcontactfile[c('mod_pos2')]))
				reciprocalust=uniprotsitetable[c(2,1)]
				colnames(reciprocalust)=colnames(uniprotsitetable)
				crosslinktable=rbind(uniprotsitetable,reciprocalust)
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
			}

			

			
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
			#old version: last value will always be match so subtract 1 to find true overlap.
			matrix.matcher=cms.xdm==cms.matrix
			matrix.matcher=matrix.matcher+0
			#matrix.overlap=sum(matrix.matcher,na.rm=TRUE)-1
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
		#
		#To allow scoring of residual score==0
		#
			max.cmscore.adj=total.crosslinks
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
		text(median(distancerange),
			max(cmscore.table.final),
			labels=paste("CMScore =",round(sum(cmscore.table.final[,2]),digits=5)))
		
		if(save.cms.table){
			question.save.table = "YES"
		} else {
			question.save.table = winDialog(type = c("yesno"),"Do you want to save the CMScore table?")
		}
		
		if(question.save.table=="YES") {
			write.table(cmscore.table.final,
				file = paste(cms.proteininput,
					"CMScoreTable",
					format(Sys.time(),"%Y%m%d%H%M%S"),
					".txt", sep=""), 
				sep="\t",
				row.names=FALSE
			)
		} else {
			print("NICE PLOT!")
		}
	}
}

#Example: cmscore(0:40,"B7I4V0")
#Once run on a single protein and range, compare to multiple structures
