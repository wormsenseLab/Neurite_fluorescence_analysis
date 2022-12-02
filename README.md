# Neurite_fluorescence_analysis
Comprehensive script for measuring fluorescence intensity in diffuse and punctate pools, inter-punctum intervals (IPI), puncta colocalization, intensity correlation quotient (ICQ) in neurites.


SAMPLE PREP: 
Grow worms at 20 degrees. 
Mount young adults on 5% agarose pads in M9.
Use 5 mM Levamisole in M9 for immobilization.

DATA COLLECTION: 
Only image the neuron closest to the cover slip. 

FILE NAMING:
Group folder name: yyyymmdd_Strain_xx \n
File prefix: yyyymmdd_Strain \n
Stitched file name: yyyymmdd_Strain_xx

IMAGE PRE-PROCESSING: 
Open image in Fiji/ImageJ. Stitched files are rotated to orient anterior end towards the RIGHT unless there are multiple worms in the frame. To distinguish between ALML/R, look for alae position relative to ALM. If alae is below TRN then ALML, otherwise ALMR. In case of ALMR, often the AVM cell body is also visible. Trace neuron using segmented line tool. Use spline fit. Always trace starting from cell body moving towards distal end. Save ROIs. Press spacebar to switch to hand mode for easy scrolling through zoomed image. Check to make sure every part of the segmented line aligns well with the neuron. Before straightening make sure to reset brightness and contrast levels. Line width for straightening: 20 px. Save straightened image. Straightened file name: yyyymmdd_Strain_xx-x Every straightened image should be easily traceable to its raw image file.

QUALITY CONTROL: Do not straighten images with the following issues (trace neuron and save ROI): \n
-- Puncta out of focus \n
-- Movement of worm during acquisition \n
-- Improper stitching

PUNCTA ANALYSIS: 
Run script Neurite_fluorescence_analysis.py.
The algorithm will ask for the following user inputs:
-- Select strain_key spreadsheet file (GUI)
-- Camera pixel width (um)
-- Select destination folder for saving analysis files (GUI)
-- Objective magnification
-- Age of worms
-- Neuron
-- Select folder where straightened images are stored (GUI)
The algorithm will create a new folder with current timestamp within the destination folder to save all output from the current run.
The algorithm will output traces of each neuron, with puncta identified as green dots. Verify visually if peaks are appropriately identified. These images are also saved as .svg files within a sufolder called "individual_traces".
The algorithm will output all the relevant dataframes as .pkl files. For all further analysis and plotting, read these .pkl in pandas and extract data as needed.
