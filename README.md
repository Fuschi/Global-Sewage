# Global-Sewage

This repository contains a network analysis of Antimicrobial Resistance Genes (ARGs) from sewage data. The analysis focuses on the identification and visualization of ARGs across various bacterial species and their interactions in sewage samples.

## Project Structure

The project is organized into the following folders:

- **analysis**: This directory contains the notebook that you can use to reproduce the analyses presented in the project.
- **data**: Here, you will find the organized database saved in R format. Due to space limitations, the raw data files are not included in this repository. The metadata on antimicrobial resistance (AMR) could be better organized and may need a thorough review.
- **plots**: All generated plots from the analyses are saved here.
- **tables**: This folder contains summary tables that summarize the network data.

## Requirements

To run the code and reproduce the analyses, you will need the **mgnet** package installed: -  (currently in development but soon to be released)  

```
library("devtools")
install_github("Fuschi/mgnet", ref = "a6642dff48cfbe5183d8f9d8bc9e983aaf5dcc05")
```
