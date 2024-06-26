# Model development for the CD16 signaling reactions



***Signaling reactions in compartments***: We consider a well-mixed system consisting of CD16 receptor, anti-CD16 antibody, CD3ζ adaptor, SFK Lck, Syk-family kinases ZAP70 and SYK, phosphatase SHP-1, and ubiquitin ligase Cbl in a simulation box representing a small region (5 μm × 5 μm × 1 μm) proximal to the cell membrane. We divided the simulation box into two compartments representing the extracellular volume (5 μm × 5 μm × 0.002 μm = 0.05 μm3) where CD16 interacts with the anti-CD16 antibody,  the plasma membrane (5×5=25 μm<sup>2</sup>), and cytosolic volume underneath the plasma membrane (5 × 5 × 1 = 25 μm<sup>3</sup>). In this model, we considered only homodimeric forms of CD3ζ consisting of 6 ITAMs.  The interactions involving the above molecules considered in the in silico model are following.

***Reactions rules***: CD16 receptors bind to the anti-CD16 antibodies to form the anti-CD16:CD16 (ligand:receptor) complexes on the plasma membrane. Then, adaptor proteins CD3ζ associate with anti-CD16:CD16 complexes at the plasma membrane. Each ITAM of CD3ζ contains two tyrosine residues that are phosphorylated by the SFK, Lck. An ITAM can be in fully unphosphorylated, or partially phosphorylated (one phosphorylated tyrosine residue), or fully phosphorylated (two phosphorylated tyrosine residues).  For simplicity we considered that each ITAM in the model can be either in an unphosphorylated (U) or fully phosphorylated (PP) state. The unphosphorylated state U in the model represents the unphosphorylated ITAM or partially phosphorylated ITAM. 

***Rules involving ZAP70***: Unbound ZAP70 molecules in the cytosol bind to the fully phosphorylated ITAMs.  ITAM-bound unphosphorylated ZAP70 is further phosphorylated by Lck kinase.  ITAM-bound ZAP70 either in phosphorylated or unphosphorylated form can unbind from the complex into the cytosol. Cbl binds to ITAM-bound ZAP70 and degrades the entire complex including CD3ζ and CD16 [^1][^2]. 

***Rules involving Syk***: Free SYK molecules (unphosphorylated or phosphorylated) in the cytosolic bind to the partially or fully phosphorylated ITAMs on CD3ζ. Because we do not have partially phosphorylated ITAMs explicitly accounted for in the model, we assumed SYK associated with the U- state of the ITAMs with a weaker affinity. ITAM-bound SYK is further phosphorylated by Lck.  Unlike ZAP70, bound unphosphorylated SYK (denoted as basally active SYK) molecules can auto-phosphorylate itself to the phosphorylated form (denoted as catalytically active SYK or pSYK) [^3] (see Fig. 2F in main text). Moreover, it can fully phosphorylate tyrosine residues on the same ITAM it is bound to or tyrosine residues on other ITAMs [^3]. Similarly, catalytically active pSYK can trans-phosphorylate ITAM bound SYK [10]. SYK also can unbind from the ITAMs to the cytosol (free). Similar to ZAP70, Cbl can associate with ITAM-bound SYK and degrade the SYK molecule along with CD16 and CD3ζ [^1].

***Rules involving SHP-1***: SHP-1 in the cytosol can bind to partially and fully phosphorylated ITAMs. Similar to the rules used for modeling SYK binding to partially phosphorylated ITAMs, we included SHP-1 binding to the U state of ITAM.  ITAM-bound SHP-1 can dephosphorylate the pZAP70/pSYK [^4][^5][^6]. SHP-1 can also unbind from the ITAMs. In our model, the abundances of SHP-1 and their affinities to ITAMs are much lesser compared to ZAP70 and SYK (see Table S1, Table S2). Thus, we did not observe a significant SHP-1 binding to ITAMs in our in silico signaling model.

***Rules involving Cbl***: Cbl is a cytokine ubiquitin ligase that binds to the ITAM-bound ZAP70 and SYK both in phosphorylated and unphosphorylated form and degrades the entire complex along with receptor CD16 and adaptor protein CD3ζ [^1][^2].

Additionally, we applied the kinetic proofreading scheme [^7] in our signaling model, which allows the breaking of all intermediate complexes into their free and unphosphorylated initial states with the rate (KPR) that is equal to the unbinding rate (koff) between the anti-CD16 antibody and CD16 receptor (see Table S1). 

***Visualization of biochemical reaction rules***: The reaction rules for rule-based model are created using the Virtual Cell software [^8]. Then V-cell code is exported in .bngl (BioNetGen) format to run in Python language using PyBioNetGen library [^9][^10]. The reaction rules can be visualized by loading .bngl file in the Virtual Cell (Vcell) software [^11]. In our model, each ITAM on CD3ζ could be in ten states. For example, ITAM can be in unphosphorylated or partially phosphorylated state (U), fully phosphorylated state (PP), ITAM(PP) bound to ZAP70 (2 states for phosphorylated and unphosphorylated state), ITAM(PP) bound to SYK (2 states for phosphorylated and unphosphorylated state), ITAM(PP) bound to SHP-1, ITAM(U) bound to SYK (phosphorylated and unphosphorylated state) and SHP-1 (3 states). This can lead to ~6<sup>10</sup> number of different reaction species which will be computationally intensive to simulate numerically using standard ordinary differential equation- or Gillespie-based applications. Therefore, we opted for a rule-based version of Gillespie, Network Free Stochastic Simulator (NFSim) to handle our extensive CD16 signaling model for WT NK cells that significantly reduces simulation runtime [27].

[^1]: Duan L, Reddi AL, Ghosh A, Dimri M, Band H. The Cbl family and other ubiquitin ligases: destructive forces in control of antigen receptor signaling. Immunity. 2004;21(1):7-17.

[^2]: Wang H-Y, Altman Y, Fang D, Elly C, Dai Y, Shao Y, et al. Cbl promotes ubiquitination of the T cell receptor ζ through an adaptor function of Zap-70. Journal of Biological Chemistry. 2001;276(28):26004-11.

[^3]: Mukherjee S, Zhu J, Zikherman J, Parameswaran R, Kadlecek TA, Wang Q, et al. Monovalent and multivalent ligation of the B cell receptor exhibit differential dependence upon Syk and Src family kinases. Science signaling. 2013;6(256).

[^4]: Das J. Activation or tolerance of natural killer cells is modulated by ligand quality in a nonmonotonic manner. Biophysical journal. 2010;99(7):2028-37.

[^5]: Mkaddem SB, Hayem G, Jönsson F, Rossato E, Boedec E, Boussetta T, et al. Shifting FcγRIIA-ITAM from activation to inhibitory configuration ameliorates arthritis. The Journal of clinical investigation. 2014;124(9):3945-59.

[^6]: Matalon O, Fried S, Ben-Shmuel A, Pauker MH, Joseph N, Keizer D, et al. Dephosphorylation of the adaptor LAT and phospholipase C–γ by SHP-1 inhibits natural killer cell cytotoxicity. Science signaling. 2016;9(429).

[^7]: McKeithan TW. Kinetic proofreading in T-cell receptor signal transduction. Proceedings of the national academy of sciences. 1995;92(11):5042-6.

[^8]: Blinov ML, Schaff JC, Vasilescu D, Moraru II, Bloom JE, Loew LM. Compartmental and spatial rule-based modeling with virtual cell. Biophysical journal. 2017;113(7):1365-72.

[^9]: Faeder JR, Blinov ML, Hlavacek WS. Rule-based modeling of biochemical systems with BioNetGen. Systems biology. 2009:113-67.

[^10]: Ali Sinan Salgam JRF. PyBioNetGen - A lightweight BioNetGen CLI
 2021 [cited 2024]. Available from: https://pybionetgen.readthedocs.io/en/latest/.

[^11]: Schaff J, Fink CC, Slepchenko B, Carson JH, Loew LM. A general computational framework for modeling cellular structure and function. Biophysical journal. 1997;73(3):1135-46.

[^12]:
