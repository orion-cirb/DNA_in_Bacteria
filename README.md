# DNA_in_bacteria

* **Developed for:** Céline
* **Team:** Espeli
* **Date:** December 2022
* **Software:** Fiji


### Images description

3D images taken with a x60 objective

3 channels:
  1. *DAPI:* DNA
  2. *phiYFP:* foci
  3. *TL phase:* bacteria

### Plugin description

* Detect bacteria on the average intensity Z-projection of channel 3 with Omnipose
* Detect DNA on the max intensity Z-projection of channel 1 with Omnipose
* In each bacterium, return distances between bacterium centroid and DNA centroid


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Omnipose** conda environment + *bact_phase_omnitorch_0* *bact_fluor_omnitorch_0* models

### Version history

Version 1 released on December 16, 2022.

