# RS_LANG_CamCAN
Repository for graph-based analyses of CamCAN resting-state data analyses

<h1 align="center">Analysis of the topological changes in the LANG Connectome across the adult lifespan</h1>

Method Overview

- Preprocessing of resting-state data using SPM12, ART
- Relying on 6mm-264-regions from Power atlas (Power et al., 2011)
- Connectivity matrices obtained via CONN Toolbox processing
- Graph-theory measures with GraphVar and BCT Toolbox on the previously defined 131 ROIs of the LANG Connectome (Roger et al., 2022)
- Hub classification, hub detection, and computation of the individual topologico-functional profile (TFP)
- Bootstrap procedure to compute the confidence intervals of log-ratios

- Investigation of the three-way interaction between AGE, TFP, & RSNs using pseudo-temporal multi-trajectory modeling (GBMTM, Magrini, 2022)
- Investigation of the most probable topological reconfiguration trajectories for each region

Highlights

- More accurate method for RSN assignment (with voxel-level precision) to ROIs 
- Definition and Justification of a Topologico-functional profile according to modular and interareal nodal centrality metrics (i.e., degree, Participation coefficient, Within-module_z-score, Betweenness centrality, Flow coefficient)
- Compositional data analysis and multivariate latent growth mixture modeling
- Shiny app to explore the most probable topological reconfiguration trajectory across AGE: https://clementguichet.shinyapps.io/shiny_camcan/
 
## Author
**Clément Guichet**
- Clement.Guichet@univ-grenoble-alpes.fr
- [ResearchGate](https://www.researchgate.net/profile/Clement-Guichet)

## Lab
- CNRS UMR 5105 - Laboratoire Psychologie et NeuroCognition (LPNC)
- https://lpnc.univ-grenoble-alpes.fr/laboratoire-0
