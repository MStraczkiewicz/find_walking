# A “one-size-fits-most” walking recognition method for smartphones, smartwatches, and wearable accelerometers.

## Purpose

We propose a walking recognition method for sub-second tri-axial accelerometer data, in which activity classification is based on the inherent features of walking: intensity, periodicity, and duration.

## Validation

x Method has been validated against 20 publicly available, annotated datasets on walking activity data collected at various body locations (thigh, waist, chest, arm, wrist).

x We demonstrated that our method can estimate walking periods with high sensitivity and specificity: average sensitivity ranged between 0.92 and 0.97 across various body locations, and average specificity for common daily activities was typically above 0.95.

x We also assessed the method’s algorithmic fairness to demographic and anthropometric variables and measurement contexts (body location, environment).

## How to use
1. Download both files (find_walking.m, find_continuous_dominant_peaks.m) to your workstation and add them to your MATLAB search path.
2. Transform your raw accelerometry data to vector magnitude.
3. Run find_walking.m as described in function description.


## Reference

Straczkiewicz, M., Huang, E.J. & Onnela, JP. A “one-size-fits-most” walking recognition method for smartphones, smartwatches, and wearable accelerometers. npj Digit. Med. 6, 29 (2023). https://doi.org/10.1038/s41746-022-00745-z
Open Access: https://rdcu.be/c6dGV 
