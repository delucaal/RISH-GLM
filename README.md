### This repository contains the code needed to perform cross-site harmonization of diffusion MRI data using the RISH-GLM method. 

In short, RISH-GLM builds on the established "*Rotational Invariant Spherical Harmonics*" ([**RISH**](https://www.sciencedirect.com/science/article/pii/S1053811918307717?via%3Dihub)) framework, and introduces the concept of cross-site harmonisation without matched training subjects. A pre-print of this method can be found [here](https://www.biorxiv.org/content/10.1101/2024.05.01.591994v2). 

### Key scripts:
- **"example_main_rish.m"**: this script creates a study-specific spatial template between two or more sites, then trains harmonization while correcting for confounding effects (e.g., age, sex). If you do not specify any covariate, this is equivalent to running the conventional RISH. This script will be maintained. 
- **example_main_rish.m**: this script creates a study-specific spatial template between two sites, and trains harmonization using the **conventional** RISH method. Please note that this script is only provided for reference and will not be maintained.  

Many thanks to my co-authors at Harvard Medical School, particularly Suheyla Cetin Karayumak, Yogesh Rathi and Ofer Pasternak for sharing their original MATLAB implementation. 

### Publications where (part) of this code was used:
- Harmonization of diffusion MRI in small vessel disease [https://pmc.ncbi.nlm.nih.gov/articles/PMC8609094/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8609094/)
- Multi-shell harmonization with RISH [https://www.sciencedirect.com/science/article/pii/S1053811922005560](https://www.sciencedirect.com/science/article/pii/S1053811922005560)


### Other resources:
- A Python implementation of the **conventional** RISH can be found here: https://github.com/pnlbwh/dMRIharmonization

