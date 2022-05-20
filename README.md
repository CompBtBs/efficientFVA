# An efficient implementation of Flux Variability Analysis for metabolic networks


![intro_paper](https://github.com/CompBtBs/efficientFVA/blob/main/Workflow.png)


Flux Variability Analysis (FVA) is an important method to analyze the range of fluxes of a metabolic network. FVA consists in performing a large number of independent optimization problems, to obtain the maximum and minimum  flux through each reaction in the network. Although several strategies to make the computation more efficient have been proposed, the computation time of an FVA can still be limiting. We present a two-step procedure to accelerate the FVA computation time that exploits the large presence within metabolic networks of sets of reactions that necessarily have an identical optimal flux value or only differ by a multiplication constant. The fist step identifies such sets of reactions. The second step computes the maximum and minimum flux value for just one element of each of set, reducing the total number of optimization problems compared to the classical FVA. 
