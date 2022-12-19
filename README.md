# reactobiome
Relative abundance of Reactions and Reactobiome from MSP table and Gene Count Table

This Python script takes as input a MSP relative abundance table (MSP against Samples)
or a gene count table, and outputs a Reaction relative abundance table and a reactobiome table.
The ID for the ouput reactions are the KEGG database reaction ids.

Usage example,
Using MSP table as input:

$ ReactionAbundance.py --gct_table  \<path to MSP table\>


Using gene count table as input:

$ ReactionAbundance.py --gct_table  \<path to input gene count table \>


The values in the reactobiome table reflect the fraction of MSP where a reaction is present 
in relation to the total number of MSP detected for that sample (times 500) and converted into int type,

Side Notes:
1) Estimation from Gene Count Table uses the Gene count table directly. The total number of MSPs
   is just calculated as the number of different MSP associated to all the genes present in a sample,
   lilkely inflating the number of MSPs.

2) Estimation from the MSP table assumes every gene associated to a MSP is present, thus if an MSP is present in the table,
   all the genes associated to it will be given the same relative abundance as the MSP, even if the gene was not present in 
   the original GeneCountTable from wich the MSP table was originally derived.
