# *_Speech Recogniser_*  
   ______
 
   C++ implementation of a speech recognition system.

   ---
   
# yesNo.cpp
 Simple code that discriminates between "YES" and "No" based on the differences in Short term energy and Zero Crossing Rate of the signal.
   
   ---
   
# vowelwriter.cpp
   
   Code to record a vowel and generate its LPC coefficients usind Levinson Durbin's algorithm and subsequently storing the cepstral coefficients.
   
   ---
   
# vowrec.cpp

Recognises a vowel based on cepstral coefficients with Tokura Distance measure. Compares against a template in the code.

---
# K_means.cpp

Implementation of K-means clustering algorithm with random initialisation. Dataset is in Universe.csv , each data point is a 12 dimensional vector of cepstral coefficients.

---
# LBG.cpp

Code for the LBG algorithm to generate codebooks of arbitrary size (multiple of 2).

---
