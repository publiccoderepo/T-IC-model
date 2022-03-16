# The "T-IC" folder contains the codes used in the paper "Temporal Cascade Model for Analyzing Spread in Evolving Networks with Disease Monitoring Applications".  


# T-IC-model  
You can generate "dtinfmax" with different interfaces separately.

# 1.Seed selection
1.1 main_seedselection.cpp
The interface for NYC, Tokyo, SG-T and Haslemere datasets to get solution sets of RSM, ESM, Greedy-IM, Max-Deg and Random. 
Example: ./dtinfmax ./data/TKY_data.csv TKY 25all_TKY 10000 25 1 20 50 10 2 1 10 0.01 10 1


1.2 main_Italy_seedselection.cpp
The interface for Italy dataset to get solution sets of RSM, ESM, Greedy-IM, Max-Deg and Random.

1.3 main_socialNetwork_seedSelection.cpp
The interface to select seed S for social network datasets.

1.4 For DIA and T-IMMU, we refer to their public codes, and the solution sets obtained from their implementation are included
in "dynamic_TIMMU_solution_set". Measuring their performance is done with:
	main_dynamic_TIMMU.cpp (NYC, Tokyo, SG-T and Haslemere datasets)
	main_BBC_dynamic_TIMMU.cpp (Haslemere datasets)
	main_Italy_dynamic_TIMMU.cpp	 (Italy datasets)


# 2. Intervention strategies
The interfaces to make spread analysis (Section 5.3)

2.1 main_spreadAnalysis.cpp (Section 5.3.1)
2.2 back_trace_nodes_BBC.cpp; back_trace_nodes_Italy.cpp; back_trace_nodes_TKYNYCSG_T.cpp (Section 5.3.2)
2.3 main_venueAnalysis.cpp (Section 5.3.3)

# 3. time_evaluate.cpp
The interface to evaluate  running efficiency

# 4. data format.
The Tokyo, NYC and SG-T datasets are organized into the format of "timeID, node_u, node_v, venueID". 
The Italy dataset is in the format of "timeID, node_u, node_v, prob".

Comment:
The output value is the original results over different measures.
Each of these results needs to be further normalized across methods.
