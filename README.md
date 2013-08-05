TF Knockdown Scripts
====================

This repository contains all the scripts required to repeat the analysis from Cusanovich et al. 2013.

Assuming that all the appropriate annotation files have been generated, here is the order of analysis scripts:

**Script Order**
1. analysis_allthreeruv2.R
2. analysis_analyzer.sh (calls analysis_KdLRTRUV2_AverageNS.R)
3. analysis_ResultsSummary.R
4. combined_combining.py
5. combined_resultsbed2factormatrix.py
6. combined_PhiTables.R
7. combined_factormatrix2results.R
8. combined_annotator.R
9. combined_venns.R
10. genome_phastcons_windows.sh (calls midpoints.py)
11. genome_greping.py
12. genome_annot_overlaps.R
13. chromatin_bed_converter.py
14. chromatin_resultsbed2factormatrix.py
15. chromatin_factormatrix2results.R