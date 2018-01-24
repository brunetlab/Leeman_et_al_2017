1. make_expressionSet_counts.R
	- reads Kallisto output (not executable)
	- annotates and collapses transcripts (executable)

2. normalization_deseq.R
	- normalizes count matrix from (1.)

3. regress_out_replicate.R
	- regresses out replicate (used for PCAs and other visualizations)

4. deseq2_testing.R
	- young vs old for each cell type, accounting for replicate