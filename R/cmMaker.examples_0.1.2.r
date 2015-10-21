#Script to create contact maps from PDB structures or from predicted structures for example files in package.
#Required file type: ".pdb" and ".txt"
#Input: cmexample("desired-distance-constraint", a pdb file and a react2 file)

map.example=function(cadistance,exampleinput) {
	if(is.numeric(cadistance)==FALSE){
		return("Error: Contact map distances must be numeric.")
	}else if(is.character(exampleinput)==FALSE){
		return("Error: Protein input must be character, add quote marks around protein name/accession.")
	}else{
			#scripttime=proc.time()		#to test runtime
			#Nullify
			x=NULL
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
			#libraries
				#library(zoo)
				#library(Rpdb)
				#library(lattice)
				#library(latticeExtra)
				#library(tcltk)
			#data
			#data(cmdata)
			#data(xlmapdata)
			
			recname=chainid=elename=NULL
			
			#import pdb from example file
			pdbfile.initial=paste("p.",exampleinput,sep="")
			pdbfile=get(pdbfile.initial)
			
			yes.pdb=c('oxa23','tbpa','ppib')
			no.pdb=c('ab57.2983','ppib.decoy')
			
			##if pdb file imported, find contacts within x angstroms and plot them to data matrix
			#first identify atoms
			xatoms=pdbfile$atoms
			
			#as.data.frame important for maintaining values in eleid and resid columns
			xmolecule=subset(as.data.frame(xatoms),recname=="ATOM")
			
			#then identify alpha-carbons in chain A
			if(any(exampleinput==yes.pdb)==TRUE) {
				xchaina=subset(xmolecule,chainid=="A")
				xalpha=subset(xchaina,elename=="CA")
			}else if(any(exampleinput==no.pdb)==TRUE){
				xalpha=subset(xmolecule,elename=="CA")
			}else {
			return("No example files with that name")
			}
				
			aamin=1
			xlstart=min(xalpha[c("resid")])
			xlend=max(xalpha[c("resid")])
			#aamax=
			
			#create matrix including any crystal structure offset in amino acid sequence
			xfill=as.data.frame(matrix(1000,nrow=(xlstart-1),ncol=ncol(xalpha)))
			colnames(xfill)=colnames(xalpha)
			xalphafilled=rbind(xfill,xalpha)
			
			#determine distance between alpha carbons (this will create the matrix of interactions)
			xdistmatrix=as.matrix(dist(xalphafilled[c("x1","x2","x3")]))
			xdistmax=max(as.matrix(dist(xalpha[c("x1","x2","x3")])))
			if(cadistance>max(xdistmax)){
				cadistance=max(xdistmax)
				#print(cadistance)				#to ensure distance input is being read correctly
				#print(max(xdistmax))			#to ensure distance max is being read correctly
			}else{
				cadistance=cadistance
				#print(cadistance)				#to ensure distance input is being read correctly
				#print(max(xdistmax))			#to ensure distance max is being read correctly
			}
			xdistmatadj=xdistmatrix[xdistmatrix>cadistance]=50
			
			###
			###
			###
			#Make crosslink contact map from example files which are 
			###
			###
			###
			
			crosslinktable=get(gsub('.decoy',"",exampleinput))
			#print(head(crosslinktable))
			crosslinkmap=as.matrix(crosslinktable[c(3,4)])
			
			###
			###
			###
			#PLOTTING contact values
			##
			##
			##
			
			colfunc1=colorRampPalette(c("lightgrey","black","darkgreen", "lightgreen"))
				colfunc2=colorRampPalette(c("white"))
				hmcols=c(colfunc1(10000*(cadistance/max(xdistmatrix))),
					colfunc2(10000*((max(xdistmatrix)-cadistance)/max(xdistmatrix)))
				)
				
			#save plot or just show it?
			question.save.plot = winDialog(type = c("yesno"),"Do you want to save the output file to PDF?")
			if(question.save.plot=="YES") {
				pdf(paste(cadistance,"angstrom_",exampleinput,".pdf",sep=""))
				print(
					levelplot(
						xdistmatrix,
						main=paste(exampleinput," Cross-link and Contact Map (",cadistance,"Angstroms)",sep=""),
						at=seq(min(xdistmatrix),max(xdistmatrix),length.out=5000),
						colorkey=list(space="right"),
						col.regions=hmcols,
						cuts=dim(array(hmcols))-1,
						ylab=paste(exampleinput," Residue Number",sep=""),
						xlab=paste(exampleinput," Residue Number",sep=""),
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
						
					) + xyplot(mod_pos1~mod_pos2|as.character(prot1),data=crosslinktable,
								col="black",
								pch=21,
								cex=1.5,
								lty=2,
								fill='magenta',
							)
				)
				dev.off()
				
			} else{
				print(levelplot(
						xdistmatrix,
						main=paste(exampleinput," Cross-link and Contact Map (",cadistance,"Angstroms)",sep=""),
						at=seq(min(xdistmatrix),max(xdistmatrix),length.out=5000),
						colorkey=list(space="right"),
						col.regions=hmcols,
						cuts=dim(array(hmcols))-1,
						ylab=paste(exampleinput," Residue Number",sep=""),
						xlab=paste(exampleinput," Residue Number",sep=""),
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
						
					) + xyplot(mod_pos1~mod_pos2|as.character(prot1),data=crosslinktable,
								col="black",
								pch=21,
								cex=1.5,
								lty=2,
								fill='magenta',
						)
				)
			}
			
			
			#print(proc.time()-scripttime)		#to test time for running each mapping
			print("NICE MAP!")
	}	
}

#cmexample(18,'oxa23')