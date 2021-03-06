# cont-3D-snGFR
Intravital microscopy in animal models is an emerging technique in life sciences with advanced applications in kidney research. Observations about morphological and functional changes of the smallest functional unit in the kidney – the nephron – can be done longitudinally throughout the progression of disease models in live  mice. In particular, it is possible to assess structural features of a nephron and directly correlate it to its glomerular fitration rate (GFR), a key parameter of kidney function.

After intravenous injection of a freely filtered, non-resorbable, fluorescent dye a time series is acquired by multiphoton laser scanning microsopy. Filtration is observed from the intraglomerular capillaries to the connected proximal tubulus and filtration is calculated after analysis of the intratubular signal intensity shift.. Translated to a simple image processing task, this could be generalized as the assessment of the flow rate in a tube. Previous methods for this analysis relied on two manually set measurement points in the tubulus and the overall tubular volume was merely estimated based on length and diameter in 2D images [1, 2]. However, the results we obtained with this approach were highly variable, especially in tubuli with high curvature.

We extended the workflow in ImageJ [3, 4] by adding continuous measurement along the entire proximal tubulus in x-y-plots of fluorescence intensity in every frame of the time series. Additionally, automatic modelling of actual tubular volume in a 3D dataset replaced the volume estimation, hence increasing robustness, accuracy and objectivity . Subsequent data analysis in R [5, 6, 7] included smoothing of x-y intensity plots, calculation of intensity shifts in every frame and normalization against tubular volume for exact assessment of single nephron GFR (snGFR) by linear regression.

[1] Kang, J. J. et al. (2006). "Quantitative imaging of basic functions in renal (patho)physiology." Am J Physiol Renal Physiol 291(2): F495-502. <br>
[2] Kidokoro, K. et al. (2019). "Evaluation of Glomerular Hemodynamic Function by Empagliflozin in Diabetic Mice Using In Vivo Imaging." Circulation 140(4): 303-315. <br>
[3] Schindelin, J. et al. (2012). "Fiji: an open-source platform for biological-image analysis." Nat Methods 9(7): 676-682. <br>
[4] Schindelin, J. et al. (2015). "The ImageJ ecosystem: An open platform for biomedical image analysis." Mol Reprod Dev 82(7-8): 518-529. <br>
[5] RStudio: Integrated Development Environment for R. R. Inc. Boston, MA, http://www.rstudio.com/. <br>
[6] R Foundation for Statistical Computing (2017). R: A Language and Environment for Statistical Computing. R Core Team. Vienna, Austria, https://www.R-project.org/. <br>
[7] Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis, Springer-Verlag New York. <br>

We acknowlegde the support of the Core Facility Cellular Imaging (CFCI) at the Medical Faculty Carl Gustav Carus, Technical University Dresden. 

The workflow was developed and tested by Friederike Kessel with the help of Hannah Kröger in the lab of Prof. Christian Hugo
Experimental Nephrology, Department of Internal Medicine of the University Hospital Carl Gustav Carus
Technical University Dresden. It is published in the NEUBIAS gateway of F1000Research:

Kessel F, Kröger H, Gerlach M et al. A new analysis approach for single nephron GFR in intravital microscopy of mice [version 1; peer review: awaiting peer review]. F1000Research 2020, 9:1372 (https://doi.org/10.12688/f1000research.26888.1)

