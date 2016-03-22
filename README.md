# Dream-olfaction-Rscript
This is a script for the DREAM Olfaction challenge SubChallenge 2.

Our sense of smell, olfaction, is the least understood sense of all and it is commonly overlooked compared to other senses. 
We rely on our sense of smell in our daily lives and it would impact our live greatly without it. 
The goal of this competition is better understand olfaction, more specific, to build a model that predict olfactory perception based on molecules’ physical and chemical properties.   
More information of the competition can be found at: 
https://www.synapse.org/#!Synapse:syn2811262/wiki/
Here I used R package Caret to buid a machine learning model that predicts the olfactory perception intensity. 

Abstract

A fundamental problem in systems neuroscience is mapping the physical properties of a stimulus to perceptual characteristics. In vision, wavelength translates into color; in audition, frequency translates into pitch. In olfaction, however, the translation of physicochemical properties of an odorant to odor perception remains unclear. The DREAM Olfaction Challenge asked competitors to build quantitative models that predict odor intensity, odor pleasantness, and 19 odor character descriptors from molecular structure using a dataset where 49 untrained human subjects rated 476 odorants. Preceding this challenge, two published models (Khan et al., 2007; Kermen et al., 2011) predicted odor pleasantness from chemical structure.  These models successfully predicted odor pleasantness in this dataset (r = 0.49, p < 0.001; r = 0.29, p < 0.001). We built models for all 21 descriptors using the Extremely Randomized Trees algorithm in R.  Our model outperformed previously published models for predicting pleasantness (r = 0.65, p < 0.001).  We found that our intensity model relied largely on physicochemical features underlying the partitioning of odors between air phase and the receptor binding pocket.  In contrast, our models of “burnt” and “bakery” relied heavily on fingerprint similarity to trained target odors as a form of template matching, rather than isolating particular physicochemical features. 
