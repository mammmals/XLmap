\name{map.example}
\alias{map.example}
\title{
Contact Map EXAMPLES: To test contact map generation based on example data.
}
\description{
Contact maps are generated from protein structures - either a model or a Protein Databank (PDB) structure) at a distances desired by the user. 
	Cross-linked site relationships for the same protein are then mapped to sites within the structure and overlaid on the structural contact map.
	The output is a pdf image of the overlaid map and a table that can be used as an input for the cmscore function.
	Example proteins:
		oxa23
		ab57.2983
		ppib
		tbpa
		ppib.decoy
	
	Example in use:
		map.example(18, oxa23)
}
\usage{
	map.example(cadistance, exampleinput)
}
\arguments{
  \item{cadistance}{
    A distance in angstroms to display the structure contact map.
}
  \item{exampleinput}{
    The name of one of the example datasets (oxa23, ab57.2983, ppib, ppib.decoy, or tbpa).
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
A PDF plot is generated with axes corresponding to the sequence of the input protein.
}
\author{
 Devin K Schweppe, Bruce Lab, University of Washington
}

\keyword{ example contact map }
\keyword{ exampleinput }

\alias{exampleinput}

