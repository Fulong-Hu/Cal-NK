# User Manual for Cal-NK.py
## Version: 1.0
## Release Date: April 18, 2025

## 1. Overview
This script provides automated calculation and analysis of optical material properties based on dielectric function (ε), reflection spectrum, and extinction coefficient. It supports deriving the refractive index (n) and extinction coefficient (k) through script automation. This document details data sources, applicability conditions, and operational procedures.

## 2. Data Sources
### 2.1 Raw Data Directories
#### 1.	epsilon(real) and epsilon(imaginary)
Generated via First-Principles Calculations, suitable for numerical analysis of the real (ε₁) and imaginary (ε₂) parts of the dielectric function
#### 2.	Reflection Spectrum
Source: Semiconductor Spectral Analysis and Fitting Calculations, Science Press (ISBN: 978-7-03-039566-5). Used for theoretical validation of reflection spectra.
#### 3.	Extinction Coefficient, Transmission Coefficient, and Multilayer Parameters
Extracted from the International Standard Refractive Index Database ([refractiveindex.info](https://refractiveindex.info/)), including extinction coefficients and multilayer optical parameters.
#### 4.	Various Spectral Scales
Example files for WS₂ material cover energy, wavelength, and wavenumber scales with both ascending/descending data orders.
Source: International Standard Refractive Index Database ([refractiveindex.info](https://refractiveindex.info/)).

## 3. Applicability and Constraints
### 3.1 Calculation Method Selection Guide
#### 1.	Dielectric Function Calculation (ε → n/k)
Use Case: Derive n and k from known ε₁ or ε₂.
Input Requirements: Ensure ε₁ (real) or ε₂ (imaginary) data is complete and properly formatted.
#### 2.	Reflectance Calculation (Reflection → n/k)
Use Case: Inversion of optical parameters for highly absorbing materials (transmittance < 4%).
Constraints: Negligible transmittance; reflection spectrum must cover the full wavelength range.
#### 3.	Transmission/Extinction Coefficient Calculation (Transmission/Extinction → n/k)
Use Case: Analysis of anti-reflective materials (reflectance ≈ 0).
Input Requirements: Extinction coefficient must align with transmission data to avoid reflectance interference.

## 4. Script Usage
### 4.1 Execution Workflow
#### 1.	Environment Setup
Copy the Cal-NK.py to the target data directory.
Install required Python libraries (NumPy, SciPy, Matplotlib, and OS).
#### 2.	Run the Script
Follow command-line prompts to select calculation mode (dielectric function/reflectance/transmission) and input parameters.
![QQ_1744975726320](https://github.com/user-attachments/assets/09c843df-a214-47e3-9e36-bbdae68d3141)
results will output as following:
![QQ_1744975871435](https://github.com/user-attachments/assets/54db3da1-2f91-4444-b7f1-7f8681ac94b3)
![image](https://github.com/user-attachments/assets/f99c1e66-f5ad-4491-86e1-d1c4a1c9302a)
#### 3.	Output Results
Generates Result.txt containing refractive index (n) and extinction coefficient (k) values.
Visualization (e.g., n/k vs. wavelength curves). Comparison Chart.jpeg

## 5. Citation Requirements
If this script is used in publications, cite as:
F. Hu, J. Qiao, X. Zhang, Y. Luo, Y. Pu, H. Hu, G. Shi, L. Li, J. Shang. Multifunctional calculator for the optical constants of refractive index and extinction coefficient

## 6. Notes
1.	Data Validation: Experimentally verify critical parameters (e.g., reflectance thresholds) to ensure compliance with physical model assumptions.
2.	Path Conventions: Avoid spaces or special characters in directory paths to prevent script errors.
3.	Version Compatibility: Requires Python 3.7.4 or higher. Upgrade older environments if needed.

Technical Support: Contact 1500148414@qq.com for assistance.

Changelog: See CHANGELOG.md.

© 2025 Cal-NK


