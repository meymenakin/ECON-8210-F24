# ECON 8210: Homework II - Quantitative Macroeconomics

Author: Mahmut Eymen Akin
Date: January 11, 2025

## Overview

This repository contains solution files for Homework II, which implements numerical methods to solve macroeconomic policy functions (consumption, labor, and capital). Methods include:
	•	Chebyshev Polynomial Projection
	•	Finite Elements Projection
	•	3rd Order Perturbation
	•	Neural Network Approximation

## Files
	1.	**econ8210hw2_meakin.pdf:** PDF report comparing the solution methods.
	2.	**econ8210hw2_meakin.tex:** LaTeX source for the report.
	3.	**econ8210hw2q1-2_meakin.ipynb:** Jupyter Notebook that contains the summary of the model, derivation of the optimality 			conditions, and Julia code that computes the policy functions
		through Chebyshev Polynomials and Finite Elements. Also saves the output.
	4.	**econ8210hw2q3_meakin.m:** MATLAB script that solves that computes policy functions through 3rd Order Perturbation. Also saves 		the output.
	5.	**econ8210hw2q3_meakin.mod:** Dynare model file for use with the perturbation MATLAB script.
	6.	**econ8210hw2q4_meakin.ipynb:** Jupyter Notebook that contains Python code that computes the policy functions through Neural 			Network Approximation.
	7.	**Homework2_2024_Fall.pdf:** PDF file containing the assignment.

## Requirements
	•	Python (for neural network)
	•	MATLAB + Dynare (for perturbation)
	•	Julia (for projection methods)

## Instructions
	•	Run .ipynb notebooks in Jupyter for Questions 1, 2, and 4.
	•	Use MATLAB + Dynare for Question 3.
	•	Compile econ8210hw2_meakin.tex to generate the report.
