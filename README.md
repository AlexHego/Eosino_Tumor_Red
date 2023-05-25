# Eosino_Tumor_Red
## Lung Cancer Tissue Analysis Script

These scripts are designed to analyze lung cancer tissue samples using the software QuPath. The script will divide the images into three regions Outer_external, Outer_internal and tumor. The script will then detect all cells and positive cells for DAB (eosinophil) in each region. The script will also export all measurements for further analysis in a .tsv file.

## Usage

To run the script, please create a project in QuPath (empty folder) and upload your images. If the image contain a tumor please create an annotation around it before to start the script.
* Download the model, unzip it and save it in the pixel classifier folder </br> 
The path must be project_folder\classifiers\pixel_classifiers\detection_tissue.json  </br>
[Download the model_for tissue detection ](https://github.com/AlexHego/Eosino_Tumor_Red/blob/main/detection_tissue.zip) </br>

* Open the script editor on QuPath. "Automate" > " Show script editor" </br>
* Open the github page of one of the scripts </br>
[ Script_Eosino_tumor ](https://github.com/AlexHego/Eosino_Tumor_Red/blob/main/Script_Eosino_Tumor_Red_chromogen.groovy) , BSD-3 licence </br>

* Copy and paste the code in the script editor
* Then, simply run the script by selecting "run" > "run for project" in QuPath.


## Output

The script will output a .tsv file for each region containing the number of cells, positive cells for DAB (eosinophil), the area and the number of positive cells per mm².

## Note

Please note that the script is designed for use with lung cancer tissue samples from the GIGA and may not be suitable for other types of tissue or analysis. It is important to validate the results with an expert before making any conclusions.



