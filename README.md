# Wideband Illumination Liquid Crystal Reconfigurable Intelligent Surfaces: Modeling, Design, and Experimental Tests
[![License: MIT + Citation Request](https://img.shields.io/badge/License-MIT--Citation-yellow.svg)](./LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2604.09214-b31b1b.svg)](https://arxiv.org/pdf/2604.09214)
[![arXiv](https://img.shields.io/badge/arXiv-2508.04331-b31b1b.svg)](https://arxiv.org/pdf/2508.04331)

This repository contains the source code for the following papers:

[1] M. Delbari, R. Neuder, A. Jiménez-Sáez, Q. Zhou, and V. Jamali, "Wideband Illumination with Liquid Crystal Reconfigurable Intelligent Surfaces: Modeling, Design, and Experimental Tests," Submitted to IEEE Transactions on Wireless Communications.

https://arxiv.org/pdf/2604.09214

[2] M. Delbari, R. Neuder, A. Jiménez-Sáez, Q. Zhou, and V. Jamali, “Near-field Liquid Crystal RIS Phase-Shift Design for Secure Wideband Illumination”, in Proc. IEEE Globecom Workshops (GC Wkshps), Taipei, Taiwan, 2025.

https://arxiv.org/pdf/2508.04331

# Usage
1. Run Fig_3.m file to generate Fig 3.
2. Run the main file, and then figures 5, 7, and 8 will be generated (Without changing parameters)
3. If you want to change parameters, you need to optimize phase shifts again. To do it:

  3.1. Download CVX from https://cvxr.com/cvx/.
  3.2. Place the CVX folder in the same folder as this project.
  3.3. Run cvx_setup.m.
  3.4. Then you can uncomment the lines in the Main.m file to optimize the RIS phase shifts.

# Acknowledgement
The work of M. Delbari, Q. Zhou, and V. Jamali was supported in part by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) within the Collaborative Research Center Multi-Mechanisms Adaptation for the Future Internet (MAKI) (SFB 1053) under Project-ID 210487104; in part by the LOEWE initiative (Hesse, Germany) within the emergenCITY Centre under Grant LOEWE/1/12/519/03/05.001(0016)/72, and in part by the German Federal Ministry for Research, Technology and Space (BMFTR) under the program of “Souverän. Digital. Vernetzt.” joint project Open6GHub plus (Project-ID 16KIS2407). Neuder and Jiménez-Sáez’s work was supported by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Project-ID 287022738 – TRR 196 MARIE within project C09. In addition, thanks goes to Merck Electronics KGaA, Darmstadt, Germany, for providing the liquid crystal mixture. 

# License and Referencing
This program is licensed under the MIT license. If you in any way use this code for research that results in publications, please cite our original article listed above. You can also use the following BibTeX entry

<pre lang="markdown">
@misc{delbari2026wideband,
  title={Wideband Illumination with Liquid Crystal Reconfigurable Intelligent Surfaces: Modeling, Design, and Experimental Tests},
  author={Delbari, Mohamadreza and Neuder, Robin and Jim{\'e}nez-S{\'a}ez, Alejandro and Zhou, Qikai and Jamali, Vahid},
  year={2026},
  url = {https://arxiv.org/pdf/2604.09214}
}
</pre>

<pre lang="markdown">
@inproceedings{delbari2025wideband,
  title={Near-field Liquid Crystal {RIS} Phase-Shift Design for Secure Wideband Illumination},
  author={Delbari, Mohamadreza and Neuder, Robin and Jim{\'e}nez-S{\'a}ez, Alejandro and Zhou, Qikai and Jamali, Vahid},
  booktitle = {Proc. IEEE Globecom},
  year={2025}
}
</pre>
