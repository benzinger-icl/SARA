# SARA

SARA: Select Atrophied Regions in Alzheimer Disease is a logistic algorithm that predicts the likelihood of a symptomatic Alzheimer disease diagnosis. It uses volumetric data determined by FreeSurfer 5.3. 

SARA additionally will plot volumetric data on Normal Aging curves and determine an individual's percentile for their age.

The results of running this script are a CSV file with all original volumetric data, though with volumes having been corrected for intracranial volume, the age-specific z-score for each subject for each region, and the SARA output: expressed both as a percentage from 0-1 and categories based upon sensitivity and specificity in a previously evaluated group of patients. If make_graphs=T, a graph of each region showing the subject's data overlayed on normal aging curves, and a graphical representation of the SARA results will additionally be created.

### Example Output Images

![Example_Normal_Aging_Curve](https://github.com/benzinger-icl/SARA/blob/master/Example_Images/Example_Right.Hippocampus_volume.png)
![Example_SARA_Image](https://github.com/benzinger-icl/SARA/blob/master/Example_Images/Example1_SARA_classification.png)

The specific data required to run SARA are: ID, Age, Intracranial volume, Left and Right Hippocampus volumes, Left and Right Amygdala volumes, Left and Right Inferior Lateral Ventricle Volume, Left and right entorinhal thickness, left and right inferior parietal thickness. 

Other thicknesses and volumes can be included if desired, but no other data (ie demographics). Region names must match those in the Example_Patient_Input.csv or the script will fail.

These are used to calculate the percent liklihood that an individual has a pattern of atrophy indicative of symptomatic Alzheimer disease. Part of this process includes adjusting all volumes (except the white matter hyperintensities if included) for intracranial volume, or head size. The information needed for this is part of the R environment that should be downloaded in addition to the R script containing the SARA_function(). 

The R environment also contains information used to plot the normal aging curves, and determine the age-specific z-score/percentile for each subject plotted on the normal aging curve. Make_graphs=True assumes the provided data is all for a single subject, and will additionally calculate the rate of change for each region over time, and the corresponding rate for the normal aging graph. 

SARA was created and tested on R version 3.5.3 and RStudio version 1.2.5019, and using data from FreeSurfer 5.3. Use of other versions may impact results or cause the script to fail.

## To Use:

	1. Download and install R and R Studio
	
	2. Download from this repository the files 
		Example_Patient_Input.csv
		SARA_Normal_Aging_Input.RData
		SARA_function.R
	
	3. Enter your FreeSurfer information into the Example_Patient_Input.csv, with each row being an MRI scan
	
	4. Open R studio and load the function by entering the following into the terminal and hitting enter
	
```{r}
source("X:/file/location/SARA_function.R")
```
	
	5. Run the function by entering the following into the terminal and hitting enter
	
```{r}
SARA(
normals_data_location ="X:/file/location/SARA_Normal_Aging_Input.RData",
patient_data_location ="X:/file/location/Example_Patient_Input.csv",
)
```
		
	6. Depending on what you want, you can add any of the following inside the parentheses, separated by commas
	
		If you want to also create the images (only use if your data is only on one person):
		
		make_graphs=TRUE
		
		If you want to only create the images and not the .csv file add this in addition to make_graphs=TRUE :
		
		make_csv=FALSE 
		
		If you want the output of the function to go somewhere other than your current working directory:
		
		output_directory = "X:/file/location/"
		
		If you want to name the columns differently than what is given in the example input, you have to specify 
		what labels you are using for the data required to run the script:
		
		ID_name = "CLINICALDATA.ID"
		
		age_label = "Age"
		
		IntraCranialVol = "IntraCranialVol"
		
		Left.Amygdala_volume = "Left.Amygdala_volume"
		
		Right.Amygdala_volume="Right.Amygdala_volume"
		
		Left.Hippocampus_volume="Left.Hippocampus_volume"
		
		Right.Hippocampus_volume="Right.Hippocampus_volume"
		
		Left.Inf.Lat.Vent_volume="Left.Inf.Lat.Vent_volume"
		
		Right.Inf.Lat.Vent_volume="Right.Inf.Lat.Vent_volume"
		
		lh_entorhinal_thickness="lh_entorhinal_thickness"
		
		rh_entorhinal_thickness="rh_entorhinal_thickness"
		
		lh_inferiorparietal_thickness="lh_inferiorparietal_thickness"
		
		rh_inferiorparietal_thickness="rh_inferiorparietal_thickness"

Creation and validation of the SARA model is outlined here: (under review, will update when paper is accepted)
SARA's Alzheimer classification has not been evaluated in those under age 45 and works best for those being evaluated for Alzheimer disease in addition to non-neurodegenerative diagnoses and/or frontotemporal dementia.
