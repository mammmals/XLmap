#adapted from read.pdb function in R-package Rpdb 2.2 to import PDBx/mmCIF files
read.pdbx.atoms=function (file) {
### To get pdbx files to work in cmmaker and scorer
    if (!file.exists(file)) 
        stop("File '", file, "'' is missing")
    lines <- readLines(file)
    recname <- substr(lines, 1, 6)
    title <- NULL
    title <- subset(lines, recname == "TITLE ")
    remark <- NULL
    remark <- subset(lines, recname == "REMARK")
    atom <- character(0)
    is.hetatm <- rep(FALSE, length(recname))
    is.atom <- rep(FALSE, length(recname))
    is.atom <- recname == "ATOM  "
    is.hetatm <- recname == "HETATM"
    atoms <- subset(lines, is.atom | is.hetatm)
	
	atoms.table=read.table(text=atoms,sep="")
	colnames(atoms.table)=c(
		'group_PDB',
		'id',
		'type_symbol',
		'label_atom_id',
		'label_alt_id',
		'label_comp_id',
		'label_asym_id',
		'label_entity_id',
		'label_seq_id',
		'pdbx_PDB_ins_code',
		'Cartn_x',
		'Cartn_y',
		'Cartn_z',
		'occupancy',
		'B_iso_or_equiv',
		'Cartn_x_esd',
		'Cartn_y_esd',
		'Cartn_z_esd',
		'occupancy_esd',
		'B_iso_or_equiv_esd',
		'pdbx_formal_charge',
		'auth_seq_id',
		'auth_comp_id',
		'auth_asym_id',
		'auth_atom_id',
		'pdbx_PDB_model_num'
	)
	return(atoms.table)
}
