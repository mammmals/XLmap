\name{cms.example}
\alias{cms.example}
\title{
CMScore EXAMPLES: To test scoring based on example data.
}
\description{
 	Contact maps are generated from example protein structures at a distances desired by the user. 
	Cross-linked site relationships for the same protein are then mapped to sites within the structure and overlaid on the structural contact map.
	At each distance, the overlap generates a residual score (lower scores [score=0] are better at each point) whereby the lowest sum of the residual scores (CMScore) 
	for multiple models corresponds to the model that best matches the observed cross-linked sites. 
	Requires '.pdb' file and text file containing cross linked sites.
	Example proteins:
		oxa23
		ab57.2983
		ppib
		tbpa
		ppib.decoy
			
	Example in use:
		cms.example(0:40, oxa23)
}
\usage{
	cms.example(distancerange, exampleinput)
}
\arguments{
  \item{distancerange}{
     A distance range in angstroms to calculate the degree of overlap over. The degree of overlap will be proportional to the distance.
}
  \item{exampleinput}{
     A protein identifier used to specify which protein to calculate a cmscore for. If only one protein exists, then input this accession. The script can grep uniprot accessions from within the "|" of protein accessions in SEQUEST outputs.
}
}
\details{
	Required inputs are:
		(1) a text file containing cross-linked protein identifications and the sites of cross-linking 
		and,
		(2) a '.pdb' file containing structural information for the queried protein.
	The text file should have at least four columns: 'prot1' 'prot2' 'mod_pos1' 'mod_pos2'
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

\keyword{ example CMScore }
\keyword{ exampleinput }

