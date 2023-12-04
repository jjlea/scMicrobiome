cd $MyDIR/scFEA-master
 
python src/scFEA.py --data_dir data \
	--input_dir $MyDIR/tables \
	--test_file hepatocyte_expr.matrix_log_filtered_1558.cells.csv \
	--moduleGene_file module_gene_m168.csv \
	--stoichiometry_matrix cmMat_c70_m168.csv \
	--sc_imputation True \
	--output_flux_file $MyDIR/flux_filtered_1558.cells.out \
	--output_balance_file $MyDIR/balance_filtered_1558.cells.out

