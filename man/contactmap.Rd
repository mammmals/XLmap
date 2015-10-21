\name{contactmap}
\alias{contactmap}
\title{
Create protein contact maps overlaid with protein cross-linking data.
}
\description{
Contact maps are generated from protein structures - either a model or a Protein Databank (PDB) structure) at a distances desired by the user. 
	Cross-linked site relationships for the same protein are then mapped to sites within the structure and overlaid on the structural contact map.
	The output is a pdf image of the overlaid map and a table that can be used as an input for the cmscore function.
}
\usage{
	contactmap(cadistance,proteininput,molfile=NULL,
	xltable=NULL,pdbx=FALSE,cm.offset=0,sites.only=FALSE,
	which.chain=NULL, save.cm.plot=TRUE, ...)
}
\arguments{
  \item{cadistance}{
    A distance in angstroms to display the structure contact map.
}
  \item{proteininput}{
    A protein identifier used to specify which protein to calculate a cmscore for. If only one protein exists, then input this accession. The script can grep uniprot accessions from within the "|" of protein accessions in SEQUEST outputs.
}
  \item{molfile}{
	A path to User's desired protein structure or model. If none is chosen the user will be prompted to find the file during the script's run.
}
  \item{xltable}{
	A path to User's desired xltable. If none is chosen the user will be prompted to find the file during the script's run.
}
  \item{pdbx}{
   Logical; Default as FALSE. To indicate whether or not the structural file being imported is a PDB, PDBx or mmCIF file.
}
  \item{cm.offset}{
   Integer specifying if an offset should be applied to XL-MS identified crosslink sites based on amino acid residue numbering within the input structure.
}
  \item{sites.only}{
   Logical; Default as FALSE. To indicate whether or not contact map will consider only the sites (sites.only=TRUE) or the protein identifiers and sites (sites.only=FALSE) from the XL-MS input.
   If your XL-MS data is from a single protein, site.only=TRUE can be used, if you need to parse out a single protein from 
}
  \item{which.chain}{
   Specifies which chain in the structure file will be used.
}
  \item{save.cm.plot}{
   Save the Contact map/XL map plot as a PDF file to the current working directory.
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
A PDF plot is generated with axes corresponding to the sequence of the input protein.
}
\references{
Schweppe et al.
}
\author{
 Devin K Schweppe, Bruce Lab, University of Washington
}

\examples{
## For mapping of Oxa23 at 18 Angstroms:
contactmap(18, 'B7I4V0', 
molfile=system.file("examples","4JF4.pdb",package="XLmap"),
xltable=system.file("examples","Oxa23.txt",package="XLmap"))
}

\keyword{ contactmap }
