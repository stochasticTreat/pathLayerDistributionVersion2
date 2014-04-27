arms = STUDY@arms

arms = loadDataArm(description="Load drug screen data",
										title="functional_drug_screen_summary", 
										scriptFile="./drug_screen_nuevo.R", 
										mainFunction=RunDrugScreen, 
										arms=arms)

arms = loadDataArm(description="Load somatic mutation data",
									 title="somatic_mutation_aberration_summary", 
									 scriptFile="somatic_mutations_processing_v8.R", 
									 mainFunction=runSomaticMutationsProcessing, 
									 arms=arms)
arms = loadDataArm(description="Load abitrary set of genes for path enrichment",
									 title="arbitrary_gene_data_input", 
									 scriptFile="generic_aberration_summary.R", 
									 mainFunction=RunGenericEnrichment, 
									 arms=arms)

arms = loadDataArm(description="Run overlap analysis",
									 title="overlap_analysis", 
									 scriptFile="OverlapAnalysisFunctions.R", 
									 mainFunction=RunOverlapAnalysis, 
									 arms=arms)

STUDY@arms = arms



