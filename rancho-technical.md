# Rancho Bioscience Technical Interview

## Goal
Create a database that provides a wholisitic report of relevant clinical information for patients participating in clinical trials by wrangling, munging, and combining data produced by various experimental methods.

## Data
- Excel file named Technical Test - Data Wrangling.xlsx

## Tools used
- python
- pandas
- openxyml
- streamlit
- jupyter

## Final products

1. Excel file with containing the final database in the sheet named "Final Database"
2. Streamlit app - app that allows the user to upload excel sheet of interest and returns the final data wrangling results (sensitive to any changes in file or sheet name)
3. Python program - program that executes the data wrangling steps and produces an updated file of the 
4. Jupyter notebook - outlines my apporach in creating the above two products

## Repo contents
1. Jupyter notebook
2. Python program
3. Streamlit app 

## Sidenote/afterthought
It might be worthwhile to differentiate between the non-numeric values in the result category. Some of them are due to very low detection levels (i.e. <LLOQ), whiles others seem to have arisen due to experimental errors (e.g. QNS). For this exercise, I followed the prompt which. instructs that differentiation should only be made between numberic and non-numeric numbers.
