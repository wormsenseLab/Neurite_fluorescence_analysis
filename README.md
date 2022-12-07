# Neurite_fluorescence_analysis
Comprehensive script for measuring fluorescence intensity in diffuse and punctate pools, inter-punctum intervals (IPI), puncta colocalization, intensity correlation quotient (ICQ) in neurites.


SAMPLE PREP:

Grow worms at 20 degrees. 

Mount young adults on 5% agarose pads in M9.

Use 5 mM Levamisole in M9 for immobilization.


DATA COLLECTION: 

Only image the neuron closest to the cover slip. 


FILE NAMING:

Group folder name: yyyymmdd_Strain_xx

File prefix: yyyymmdd_Strain

Stitched file name: yyyymmdd_Strain_xx


IMAGE PRE-PROCESSING: 

Open image in Fiji/ImageJ. Stitched files are rotated to orient anterior end towards the RIGHT unless there are multiple worms in the frame. To distinguish between ALML/R, look for alae position relative to ALM. If alae is below TRN then ALML, otherwise ALMR. In case of ALMR, often the AVM cell body is also visible. Trace neuron using segmented line tool. Use spline fit. Always trace starting from cell body moving towards distal end. Save ROIs. Press spacebar to switch to hand mode for easy scrolling through zoomed image. Check to make sure every part of the segmented line aligns well with the neuron. Before straightening make sure to reset brightness and contrast levels. Line width for straightening: 20 px. Save straightened image. Straightened file name: yyyymmdd_Strain_xx-x Every straightened image should be easily traceable to its raw image file.


QUALITY CONTROL: 

Do not straighten images with the following issues:

-- Puncta out of focus

-- Movement of worm during acquisition

-- Improper stitching


PUNCTA ANALYSIS: 

Download all the files and directories in the repository.

Use the environment.yml file to create a new environment for running this script. 

  For a quick guide on how to create a new environment from a .yml file, go to:

  https://docs.anaconda.com/navigator/tutorials/manage-environments/#importing-an-environment (for Anaconda Navigator)

  https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file (for Miniconda)

Modify the strain_key.xlsx file to add details of the strain(s) to be analyzed (if it is not already present) and save it. The spreadsheet is pre-populated with the details of all the strains imaged and analyzed in the manuscript. To add a new entry, please follow the structure and syntax for each field closely to avoid errors. 

Run script Neurite_fluorescence_analysis.py.

The algorithm will ask for the following user inputs:

-- Select strain_key spreadsheet file (Select the "strain_key.xlsx" file in the pop-up File dialog box)

-- Select destination folder for saving analysis files (Select the "sample_output" directory in the pop-up File dialog box for a demo run)

-- Select folder where straightened images are stored (Select the "sample_data" directory in the pop-up File dialog box for demo output)

-- Camera pixel width (Enter the value in microns in the console window: for sample data the value is 7.54)

-- Objective magnification (Enter the value used in the console window: 60x for sample data)

-- Age of worms (Enter the age description in the console window: Adult for sample data)

-- Neuron (Enter the neuron name in the console window: ALM for sample data)

The algorithm will create a new folder with current timestamp within the destination folder to save all output from the current run.

The algorithm will output traces of each neuron, with puncta identified as green dots. Verify visually if peaks are appropriately identified. These images are also saved as .png files within a subdiredctory called "individual_traces".

The algorithm will output all the relevant dataframes as .pkl files. For all further analysis and plotting, read these .pkl in pandas and extract data as needed.

For more details about the code read the Wiki pages linked to this repository.
