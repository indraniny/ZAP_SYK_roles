# Model development for the CD16 signaling reactions



***Signaling reactions in compartments***: We consider a well-mixed system consisting of CD16 receptor, anti-CD16 antibody, CD3ζ adaptor, SFK Lck, Syk-family kinases ZAP70 and SYK, phosphatase SHP-1, and ubiquitin ligase Cbl in a simulation box representing a small region (5 μm × 5 μm × 1 μm) proximal to the cell membrane. We divided the simulation box into two compartments representing the extracellular volume (5 μm × 5 μm × 0.002 μm = 0.05 μm3) where CD16 interacts with the anti-CD16 antibody,  the plasma membrane (5×5=25 μm<sup>2</sup>), and cytosolic volume underneath the plasma membrane (5 × 5 × 1 = 25 μm<sup>3</sup>). In this model, we considered only homodimeric forms of CD3ζ consisting of 6 ITAMs.  The interactions involving the above molecules considered in the in silico model are following.

***Reactions rules***: CD16 receptors bind to the anti-CD16 antibodies to form the anti-CD16:CD16 (ligand:receptor) complexes on the plasma membrane. Then, adaptor proteins CD3ζ associate with anti-CD16:CD16 complexes at the plasma membrane. Each ITAM of CD3ζ contains two tyrosine residues that are phosphorylated by the SFK, Lck. An ITAM can be in fully unphosphorylated, or partially phosphorylated (one phosphorylated tyrosine residue), or fully phosphorylated (two phosphorylated tyrosine residues).  For simplicity we considered that each ITAM in the model can be either in an unphosphorylated (U) or fully phosphorylated (PP) state. The unphosphorylated state U in the model represents the unphosphorylated ITAM or partially phosphorylated ITAM. 

***Rules involving ZAP70***: Unbound ZAP70 molecules in the cytosol bind to the fully phosphorylated ITAMs.  ITAM-bound unphosphorylated ZAP70 is further phosphorylated by Lck kinase.  ITAM-bound ZAP70 either in phosphorylated or unphosphorylated form can unbind from the complex into the cytosol. Cbl binds to ITAM-bound ZAP70 and degrades the entire complex including CD3ζ and CD16 [18, 19]. 

***Rules involving Syk***: Free SYK molecules (unphosphorylated or phosphorylated) in the cytosolic bind to the partially or fully phosphorylated ITAMs on CD3ζ. Because we do not have partially phosphorylated ITAMs explicitly accounted for in the model, we assumed SYK associated with the U- state of the ITAMs with a weaker affinity. ITAM-bound SYK is further phosphorylated by Lck.  Unlike ZAP70, bound unphosphorylated SYK (denoted as basally active SYK) molecules can auto-phosphorylate itself to the phosphorylated form (denoted as catalytically active SYK or pSYK) [10] (see Fig. 2F in main text). Moreover, it can fully phosphorylate tyrosine residues on the same ITAM it is bound to or tyrosine residues on other ITAMs [10]. Similarly, catalytically active pSYK can trans-phosphorylate ITAM bound SYK [10]. SYK also can unbind from the ITAMs to the cytosol (free). Similar to ZAP70, Cbl can associate with ITAM-bound SYK and degrade the SYK molecule along with CD16 and CD3ζ [18].

***Rules involving SHP-1***: SHP-1 in the cytosol can bind to partially and fully phosphorylated ITAMs. Similar to the rules used for modeling SYK binding to partially phosphorylated ITAMs, we included SHP-1 binding to the U state of ITAM.  ITAM-bound SHP-1 can dephosphorylate the pZAP70/pSYK [14, 20, 21]. SHP-1 can also unbind from the ITAMs. In our model, the abundances of SHP-1 and their affinities to ITAMs are much lesser compared to ZAP70 and SYK (see Table S1, Table S2). Thus, we did not observe a significant SHP-1 binding to ITAMs in our in silico signaling model.

***Rules involving Cbl***: Cbl is a cytokine ubiquitin ligase that binds to the ITAM-bound ZAP70 and SYK both in phosphorylated and unphosphorylated form and degrades the entire complex along with receptor CD16 and adaptor protein CD3ζ [18, 19] .

Additionally, we applied the kinetic proofreading scheme [22] in our signaling model, which allows the breaking of all intermediate complexes into their free and unphosphorylated initial states with the rate (KPR) that is equal to the unbinding rate (koff) between the anti-CD16 antibody and CD16 receptor (see Table S1). 

***Visualization of biochemical reaction rules***: The reaction rules for rule-based model are created using the Virtual Cell software [23]. Then V-cell code is exported in .bngl (BioNetGen) format to run in Python language using PyBioNetGen library [24, 25]. The reaction rules can be visualized by loading .bngl file in the Virtual Cell (Vcell) software [26]. In our model, each ITAM on CD3ζ could be in ten states. For example, ITAM can be in unphosphorylated or partially phosphorylated state (U), fully phosphorylated state (PP), ITAM(PP) bound to ZAP70 (2 states for phosphorylated and unphosphorylated state), ITAM(PP) bound to SYK (2 states for phosphorylated and unphosphorylated state), ITAM(PP) bound to SHP-1, ITAM(U) bound to SYK (phosphorylated and unphosphorylated state) and SHP-1 (3 states). This can lead to ~610 number of different reaction species which will be computationally intensive to simulate numerically using standard ordinary differential equation- or Gillespie-based applications. Therefore, we opted for a rule-based version of Gillespie, Network Free Stochastic Simulator (NFSim) to handle our extensive CD16 signaling model for WT NK cells that significantly reduces simulation runtime [27].

