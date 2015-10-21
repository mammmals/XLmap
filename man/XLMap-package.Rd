\name{XLMap-package}
\alias{XLMap-package}
\alias{XLMap}
\docType{package}
\title{
Map protein cross-links sites to protein structure contact maps.
}
\description{
Package designed to overlay a protein structure contact map and residues cross-linked within a protein
	for visualization of cross-linked residues and to calculate the agreement between these two maps.
   to quantify the degree of overlap.
}
\details{
\tabular{ll}{
	Package: \tab XLMap\cr
	Type: \tab Package\cr
	Version: \tab 0.1.4\cr
	Date: \tab 2015-08-12\cr
	License: \tab GPL-2\cr
}
	Required inputs are:
		(1) a text file containing cross-linked protein identifications and the sites of cross-linking 
		and,
		(2) a '.pdb' file containing structural information for the queried protein.
	The text file should have at least four columns: 'prot1' 'prot2' 'mod_pos1' 'mod_pos2'
	These specify the protein to be used (e.g. by Uniprot Accession, 'P23284') and the
	positions within that protein that are interacting. As long as the columns are specified
	as described, contact map generation and scores can parse the input protein accession
	from tables containing interactions within and amongst multiple proteins.
}
\author{
Devin K Schweppe, PhD \cr
Maintainer: Devin K Schweppe <dkschwep@uw.edu>
}
\references{
Schweppe et al.
}
\keyword{ contactmap }
\keyword{ cmscore }
\keyword{ protein cross-linking }
