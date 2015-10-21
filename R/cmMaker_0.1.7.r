contactmap=function(cadistance,proteininput, molfile=NULL,xltable=NULL, 
					pdbx=FALSE,cm.offset=0,sites.only=FALSE,which.chain=NULL, 
					save.cm.plot=TRUE, ...) {
	#scripttime=proc.time()		#to test runtime
	if(is.numeric(cadistance)==FALSE){
		return("Error: Contact map distances must be numeric.")
	}else if(is.character(proteininput)==FALSE){
		return("Error: Protein input must be character, add quote marks around protein name/accession.")
	}else{
		#Nullify
		x=NULL
		pdbfile=NULL
		xatoms=NULL
		xmolecule=NULL
		pdbx.atoms=NULL
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
		xdistmax=NULL
		xfill=NULL
		xalphafilled=NULL
		xdistmatadj=NULL
		
		recname=chainid=elename=NULL
		
		if(is.null(which.chain)) {
			which.chain="A"
		} else {
			which.chain=which.chain
		}
		
		#libraries
			#library(zoo)
			#library(Rpdb)
			#library(lattice)
			#library(latticeExtra)
			#library(tcltk)

		#import pdb or pdbx file
		if(pdbx==FALSE) {
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
			colnames(pdbx.atoms)=c(
				'recname',#'group_PDB',
				'id',
				'type_symbol',
				'elename',#'label_atom_id',
				'label_alt_id',
				'label_comp_id',
				'chainid',#'label_asym_id',
				'label_entity_id',
				'label_seq_id',
				'pdbx_PDB_ins_code',
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
				'pdbx_formal_charge',
				'resid',#'auth_seq_id',
				'auth_comp_id',
				'auth_asym_id',
				'auth_atom_id',
				'pdbx_PDB_model_num'
			)
			xmolecule=pdbx.atoms
		}
		
		
		#then identify alpha-carbons in the user defined chain
		xchaina=subset(xmolecule,chainid==which.chain)
		xalpha=subset(xchaina,elename=="CA")
			
		aamin=1
		xalpha.resid=xalpha[c("resid")]
		xalpha.resid[]=lapply(xalpha.resid,as.numeric)
		xlstart=min(xalpha.resid)
		xlend=max(xalpha.resid)
		#aamax=
		
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
			#print(cadistance)			#to ensure distance input is being read correctly
			#print(max(xdistmax))		#to ensure distance max is being read correctly
		}else{
			cadistance=cadistance
			#print(cadistance)			#to ensure distance input is being read correctly
			#print(max(xdistmax))		#to ensure distance max is being read correctly
		}
		xdistmatadj=xdistmatrix[xdistmatrix>cadistance]=50
		
		###
		###
		###
		#Make crosslink contact map from react2 table
		###
		###
		###
		
		#input is uniprot ID
		#print("Choose React2 File for Crosslink Sites")
		if(is.null(xltable)){
			xlcontactfile=read.delim(
				tk_choose.files(caption="Choose text file containing crosslinked relationships."),as.is=TRUE)
		} else {
			xlcontactfile=read.delim(xltable,as.is=TRUE)
		}
		
		#input as uniprot accession or as sequest output
		if(sites.only==FALSE) {
			if(grepl("|",xlcontactfile[1,c('prot1')],fixed=TRUE)){
		#Find Uniprot Accession number between pipes and and site of modification then output the results
				upaccs1=sapply(xlcontactfile[c('prot1')],function(y) {
							substring(
								y,
								first = 4,
								last = gregexpr("\\Q|\\E.",y)[[1]][c(2)]-1
							)
					}
				)

		#Find Uniprot Accession number between pipes and and site of modification then output the results
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
		}
		
		#Creates intraprotein link table based on proteininput
		#
		# AND based on cm.offset
		#
		xlcontactfile[c('mod_pos1')]=xlcontactfile[c('mod_pos1')]+cm.offset
		xlcontactfile[c('mod_pos2')]=xlcontactfile[c('mod_pos2')]+cm.offset
		if(sites.only==FALSE) {
			uniprotsitetable=as.data.frame(
				cbind(upaccs1,upaccs2,xlcontactfile[c('mod_pos1')],xlcontactfile[c('mod_pos2')])
			)
			#Double the subsetting
			#print(proteininput)
			upinputselect1=subset(uniprotsitetable,upaccs1==proteininput)
			upinputselect2=subset(upinputselect1,upaccs2==proteininput)
			reciprocalust=upinputselect2[c(1,2,4,3)]
			colnames(reciprocalust)=colnames(upinputselect2)
			crosslinktable=rbind(upinputselect2,reciprocalust)
		} else {
			uniprotsitetable=as.data.frame(cbind(xlcontactfile[c('mod_pos1')],xlcontactfile[c('mod_pos2')]))
			reciprocalust=uniprotsitetable[c(2,1)]
			colnames(reciprocalust)=colnames(uniprotsitetable)
			crosslinktable=rbind(uniprotsitetable,reciprocalust)
		}
		
		###
		###
		###
		#PLOTTING contact values
		##
		##
		##
		
		colfunc1=colorRampPalette(c("lightgrey","black","darkgreen", "lightgreen"))
			colfunc2=colorRampPalette(c("white"))
			hmcols=c(colfunc1(10000*(cadistance/max(xdistmatrix,na.rm=TRUE))),
				colfunc2(10000*((max(xdistmatrix,na.rm=TRUE)-cadistance)/max(xdistmatrix,na.rm=TRUE)))
			)
			
		#save plot or just show it?
		if(save.cm.plot){
			question.save.plot = "YES"
		} else {
			question.save.plot = winDialog(type = c("yesno"),"Do you want to save the output file to PDF?")
		}
		
		if(question.save.plot=="YES") {
			pdf(paste(cadistance,"angstrom_",proteininput,".pdf",sep=""))
			print(
				levelplot(
					xdistmatrix,
					main=paste(proteininput," Cross-link and Contact Map (",cadistance,"Angstroms)",sep=""),
					at=seq(min(xdistmatrix,na.rm=TRUE),max(xdistmatrix,na.rm=TRUE),length.out=5000),
					colorkey=list(space="right"),
					col.regions=hmcols,
					cuts=dim(array(hmcols))-1,
					ylab=paste(proteininput," Residue Number",sep=""),
					xlab=paste(proteininput," Residue Number",sep=""),
					#make axes legible
					scales=list(
						y=list(
							limits=c(1,xlend),
							at=pretty(seq_len(ncol(xdistmatrix)))
						),
						x=list(
							limits=c(1,xlend),
							at=pretty(seq_len(ncol(xdistmatrix)))
						)
					)
					#Use ContactMapMaker to overlay crosslinked sites with PDB contact map
					
				) + xyplot(mod_pos1~mod_pos2,data=crosslinktable,
							col="black",
							pch=21,
							cex=1.5,
							lty=2,
							fill='magenta',
						)
						
			)
			dev.off()
			write.table(crosslinktable,
				file = paste(proteininput,"CrosslinkTable",format(Sys.time(),"%Y%m%d%H%M%S"),".txt", sep=""),
				sep="\t",
				row.names=F
			)
			
		} else{
			print(levelplot(
					xdistmatrix,
					main=paste(proteininput," Cross-link and Contact Map (",cadistance,"Angstroms)",sep=""),
					at=seq(min(xdistmatrix,na.rm=TRUE),max(xdistmatrix,na.rm=TRUE),length.out=5000),
					colorkey=list(space="right"),
					col.regions=hmcols,
					cuts=dim(array(hmcols))-1,
					ylab=paste(proteininput," Residue Number",sep=""),
					xlab=paste(proteininput," Residue Number",sep=""),
					#make axes legible
					scales=list(
						y=list(

							limits=c(1,xlend),
							at=pretty(seq_len(ncol(xdistmatrix)))
						),
						x=list(
							limits=c(1,xlend),
							at=pretty(seq_len(ncol(xdistmatrix)))
						)
					)
					#Use ContactMapMaker to overlay crosslinked sites with PDB contact map
					
				) + xyplot(mod_pos1~mod_pos2,data=crosslinktable,
							col="black",
							pch=21,
							cex=1.5,
							lty=2,
							fill='magenta',
					)
			)
			write.table(crosslinktable,
				file = paste(proteininput,"CrosslinkTable",format(Sys.time(),"%Y%m%d%H%M%S"),".txt", sep=""),
				sep="\t",
				row.names=F
			)
			}
		
		#print(proc.time()-scripttime) 			#to test time for running each mapping
		print("NICE MAP!")
	}
}

#contactmap(cadistance=18,proteininput='B7I4V0',pdbx=FALSE,cm.offset=0,sites.only=FALSE,which.chain="A")