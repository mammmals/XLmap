\name{cmscore}
\alias{cmscore}
\title{
Score contact maps based on intersection with protein contacts.
}
\description{
 	Contact maps are generated from protein structures - either a model or a Protein Databank (PDB) structure) at a distances desired by the user. 
	Cross-linked site relationships for the same protein are then mapped to sites within the structure and overlaid on the structural contact map.
	At each distance, the overlap generates a residual score (lower scores [score=0] are better at each point) whereby the lowest sum of the residual scores (CMScore) 
	for multiple models corresponds to the model that best matches the observed cross-linked sites. 
	Requires '.pdb' file and text file containing cross linked sites.
}
\usage{
	cmscore(distancerange, cms.proteininput,molfile=NULL,
				xltable=NULL,cms.pdbx=FALSE,cms.offset=0,
				cms.sites.only=FALSE,cms.which.chain=NULL,
				save.cms.table=TRUE, ...)
}
\arguments{
  \item{distancerange}{
     A distance range in angstroms to calculate the degree of overlap over. The degree of overlap will be proportional to the distance.
}
  \item{cms.proteininput}{
     A protein identifier used to specify which protein to calculate a cmscore for. If only one protein exists, then input this accession. The script can grep uniprot accessions from within the "|" of protein accessions in SEQUEST outputs.
}
\item{molfile}{
	A path to User's desired protein structure or model. If none is chosen the user will be prompted to find the file during the script's run.
}
\item{xltable}{
	A path to User's desired xltable. If none is chosen the user will be prompted to find the file during the script's run.
}
  \item{cms.pdbx}{
   Logical; Default as FALSE. To indicate whether or not the structural file being imported is a PDB, PDBx or mmCIF file.
}
  \item{cms.offset}{
   Integer specifying if an offset should be applied to XL-MS identified crosslink sites based on amino acid residue numbering within the input structure.
}
  \item{cms.sites.only}{
   Logical; Default as FALSE. To indicate whether or not contact map will consider only the sites (sites.only=TRUE) or the protein identifiers and sites (sites.only=FALSE) from the XL-MS input.
   If your XL-MS data is from a single protein, site.only=TRUE can be used, if you need to parse out a single protein from 
}
  \item{cms.which.chain}{
   Specifies which chain in the structure file will be used.
}
  \item{save.cms.table}{
   Save the CMScore table used to generate R plots and the final CMScore to the current working directory.
}
  \item{...}{
   Additional arguments.
}
}

\details{
	Required inputs are:
		(1) a text file containing cross-linked protein identifications (optional) and the sites of cross-linking 
		and,
		(2) a PDB, PDBx or mmCIF file containing structural information for the queried protein.
	
	The text file should have these column names: 
		site.only=FALSE:'prot1' 'prot2' 'mod_pos1' 'mod_pos2'
		site.only=TRUE:'mod_pos1' 'mod_pos2'
	These specify the protein to be used (e.g. by Uniprot Accession, 'P23284') and the positions within that protein that are interacting.
	As long as the columns are specified as described, contact map generation and scores can parse the input protein accession from tables containing
	interactions within and amongst multiple proteins.
	The protein input needs to be quoted to ensure that it is input as a character.
}
\value{
	Returns a plot and the CMScore overlaid for the desired protein and distance range.
}
\references{
 Schweppe et al. 
}
\author{
	Devin Schweppe, Bruce Lab, University of Washington
}

\examples{
## For scoring of Oxa23 at 1:40 Angstroms:
cmscore(1:40, 'B7I4V0', 
	molfile=system.file("examples","4JF4.pdb",package="XLmap"),
	xltable=system.file("examples","Oxa23.txt",package="XLmap"))
}
\keyword{ cmscore }
\keyword{ residual score }
