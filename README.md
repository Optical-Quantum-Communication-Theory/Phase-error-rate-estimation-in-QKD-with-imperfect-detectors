# Phase error rate estimation in QKD with imperfect detectors

This is a public version of the code used in *Phase error rate estimation in QKD with imperfect detectors* \[[arXiv link](https://arxiv.org/abs/2408.17349)\]. This was built for [v2.0.2](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.0.2) of the Open QKD Security package.

This code plots key rates for the decoy-state BB84 protocol in the presence of detector imperfections such as unequal detection efficiencies and dark count rates.


## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

### As zip
1. Download the linked version of the code from above and follow all [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/v2.0.2).
2. Also follow the additional Mosek install instructions if you want an exact match.
3. Download the latest release on the side bar and unzip in your preferred directory and add this folder to the Matlab path.
4. Run mainFigOne.m, maainFigTwo.m etc for the relevant plots.


### with git
1. Clone this repository and its exact submodules navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/Phase-error-rate-estimation-in-QKD-with-imperfect-detectors
```
2. Follow all further [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/v2.0.2).
3. Also follow the additional Mosek install instructions if you want an exact match.
3. Add this folder to the MATLAB path
4. Run mainFigOne.m, maainFigTwo.m etc for the relevant plots.
