# DNA_in_bacteria

* **Developed for:** Céline
* **Team:** Espeli
* **Date:** January 2023
* **Software:** Fiji


### Images description

3D images taken with a x60 objective

2 channels:
  1. *DAPI:* DNA
  2. *TL phase:* bacteria

### Plugin description

* Detect bacteria on the average intensity Z-projection of channel 2 with Omnipose
* Detect DNA on the average intensity Z-projection of channel 1 with Omnipose
* In each bacterium, return distances between bacterium centroid and DNA centroid


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Omnipose** conda environment + *bact_phase_omnitorch_0* and *bact_fluor_omnitorch_0* models

### Version history

Version 1 released on January 5, 2023.
Version 2 released on July 4, 2023: plugin able to handle images with multiple series and times.

