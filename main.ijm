
run("Close All");
setBatchMode(false);


// Declare global variable
var NR2F1_image; // NR2F1 image name
var voronoi_roi_save_path;
var DAPI_image_pathvar; 
var iter = 5; 
// NUCLEI MASK

// NUCLEI MASK Create
DAPI_image_path = "D:/Triple-Double-Cell-Detect/data/1-DAPI.jpg" ; // NR2F1 image path for stroma
nucleus_mask = create_nuclei_mask(DAPI_image_path); //function:: CREATE THE NUCLEUS MASK
// 001_Nucleus_Mask

// ---------------------------------------------------
// NUCLEI MASK FROM DAPI FUNCTION
function create_nuclei_mask(DAPI_image_path) {
// FUNCTION TO CREATE A NUCLEUS MASK
// INPUT: DAPI image
// OUTPUT: NUCLEUS Mask image title
// e.g. INPUT: 20A_Tu01_Field01_DAPI.tif
// e.g. OUTPUT: 001 Nucleus Mask

open(DAPI_image_path);
DAPI_info = getTitle();
selectImage(DAPI_info);

// Create a temporary folder for ROI
path = getDirectory("current");
// Voronoi roi save path
voronoi_roi_save_path = path + 'roi';
// Check if the folder exists; if not, create it
if (!File.exists(voronoi_roi_save_path)) {
File.makeDirectory(voronoi_roi_save_path);
print("Created output folder: " + voronoi_roi_save_path);
} else {
print("Output folder already exists: " + voronoi_roi_save_path);
}

// Get the Nucleus(DAPI) fluorescent channel image
selectWindow(DAPI_info);

run("Split Channels");
green_DAPI_info = DAPI_info + ' (green)';
blue_DAPI_info = DAPI_info + ' (blue)';
red_DAPI_info = DAPI_info + ' (red)';

// CLOSE THE GREEN AND RED CHANNEL
selectImage(green_DAPI_info);
close();
selectImage(red_DAPI_info);
close();

// SELECT THE BLUE CHANNEL
selectImage(blue_DAPI_info);
run("Blue"); // Change LUT to blue
setMinAndMax(0, 35); // Set min and max for contrast
rename("DAPI_Fluorescent_intensity_image");

// StarDist
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'DAPI_Fluorescent_intensity_image', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.479071', 'nmsThresh':'0.0', 'outputType':'Both', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'true', 'showProbAndDist':'false'], process=[false]");

// ZOOM OUT LABEL IMAGE
selectImage("Label Image");
run("Out [-]");
run("Out [-]");

// SUPERIMPOSE THE ROI TO FLUORESCENT IMAGE
selectImage('DAPI_Fluorescent_intensity_image');
run("From ROI Manager"); // Get the overlay in the image

// CREATE MASK
run("Create Mask");
run("Watershed");
rename("001_Nucleus_Mask");

//SAVE THE NUCLEUS ROIS
roiManager("reset");
close("ROI Manager");

run("Out [-]");
run("Out [-]");

// close images
close("DAPI_Fluorescent_intensity_image"); // Fluorescent intensity image
close("Label Image"); // Label image
return getTitle();
}
// ---------------------------------------------------

// CELL MASK USING VORONOI
sox9_mask_info = create_voronoi_sox9_mask(nucleus_mask);

// ---------------------------------------------------
function create_voronoi_sox9_mask(nucleus_mask_info){
// FUNCTION TO CREATE A cell mask using voronoi method without subtraction
// OPERATION: Output a label cells, with area > 150.
// INPUT: Nucleus mask image info
// OUTPUT: Cell Mask image title 004_voronoi_mask
// e.g. INPUT: 001_Nucleus_Mask
// e.g. OUTPUT: 004_voronoi_Cell_Mask

// SELECT NUCLEUS MASK
selectWindow(nucleus_mask_info);
run("Duplicate...", "title=" + "temp_Cell_Mask"); // Duplicate the cell mask

// RUN VORONOI DIAGRAM
selectWindow("temp_Cell_Mask");
run("Voronoi");
run("Enhance Contrast...", "saturated=0.3 normalize equalize");
setAutoThreshold("Default dark");
setThreshold(0, 0);
run("Convert to Mask");
rename("004_voronoi_Cell_Mask");

// Zoom out
run("Out [-]");
run("Out [-]");

// Analyze Particle Filter (150 - Infinity) to consider for cell
selectWindow("004_voronoi_Cell_Mask");
roiManager("reset");
run("Analyze Particles...", "size=150-Infinity include add");
roiManager("Deselect");

// SAVE THE ROIs
selectWindow("004_voronoi_Cell_Mask"); // Select original Voronoi image

voronoi_roi_save_path = voronoi_roi_save_path + "/voronoiRoi.zip";
// voronoi_roi_save_path = "D:/Maddie/01-18-2024/Debug/Debug_4/Output/voronoi.zip";
roiManager("save", voronoi_roi_save_path); // save the ROIs

run("Remove Overlay"); // Remove Overlay
return getTitle();
}
// ---------------------------------------------------

// Write the function here 


// ---------------------------------------------------
// LOAD ROIS AND ASSIGN EMPTY ARRAYS

roiManager("open", voronoi_roi_save_path); // Load ROIs
run("Clear Results"); // Clear results
total_cells_counts = roiManager("count"); // Total count of cells in the ROIs

// Define enpty arry to store the values for cell and nuclei measurements
nr2f1_nuclei_mean_intensity = newArray(total_cells_counts); // nr2f1 nuclei mean intensity
nr2f1_nuclei_area = newArray(total_cells_counts); // nr2f1 nuclei area
nr2d1_nuclei_x = newArray(total_cells_counts); // nr2f1 nuclei x coordinate
nr2f1_nuclei_y = newArray(total_cells_counts); // nr2f1 nuclei  y cooreindate

sox9_mean_intensity = newArray(total_cells_counts); // cell mean intensity
sox9_area = newArray(total_cells_counts); // Cell area
sox9_x = newArray(total_cells_counts); // cell x coordinate
sox9_y = newArray(total_cells_counts); // cell y coordinate

// LOAD SOX9 AND NR2F1 RAW IMAGE

// sox9 --> Entire Cell 
image_path = "D:/Triple-Double-Cell-Detect/data/3-SOX9.jpg" ; // NR2F1 image path
open(image_path); 
rename("SOX9"); // rename the image 

// nr2f1 --> nuclei only
image_path = "D:/Triple-Double-Cell-Detect/data/2-NR2F1.jpg" ; // NR2F1 image path
open(image_path); 
rename("NR2F1"); // rename the image 

roiManager("reset"); // Reset ROI manager 
// For loop to go over all the cells detected with Voronoi method
//for (i = 0; i < total_cells_counts; i++) {
for (i = 0; i < iter; i++) {
wait(1000);
print("Iteration...", i+1);// loop iteration
roiManager("open", voronoi_roi_save_path); // Load the ROIs original 


// SOX 9 MEASUREMENTS -- > Entire cell 

selectWindow("SOX9"); // select the voronoi image
roiManager("Select", i); // Select single cell in the ROI
roiManager("measure"); // Measure
sox9_area_val = getResult('Area', 0); // Area of a single cell
sox9_mean_intensity_val = getResult('Mean', 0); // Mean Intensity of a single cell
sox9_X_val = getResult('X', 0); // X cell measurement
sox9_Y_val = getResult('Y', 0); // Y cell measurement

// STORE THE RESULT for Cell measurement
sox9_area[i] = sox9_area_val; // Cell area
sox9_mean_intensity[i] = sox9_mean_intensity_val; // cell mean intensity
sox9_x[i] = sox9_X_val; // cell x coordinate
sox9_y[i] = sox9_Y_val; // cell y coordinate

// NR2F1 NUCLEAR MEASUREMENT 
wait(1000);
// Duplicate it 
selectImage("001_Nucleus_Mask"); // Select the Image
run("Duplicate...", "title=" + "Nucleus"); // Duplicate the image

// Select the Nuclei mask
selectWindow("Nucleus"); // select the voronoi image

// Create a single mask from the ROI manager
roiManager("Select", i); // Select single cell in the ROI
selectWindow("Nucleus"); // select the voronoi image
run("Clear Outside"); // Clear outside
run("Clear Results"); // clear the results
run("Create Mask"); // Create a mask
roiManager("reset"); // Reset ROI


selectWindow("Nucleus");
run("Analyze Particles...", "size=150-Infinity include add");
//roiManager("Deselect");
close("Nucleus");
run("Clear Results"); // Clear results 

// Overlay the ROI on original NR2F1 image to Measure the intensity

selectWindow("NR2F1"); // select the voronoi image
roiManager("Deselect");

// Ques. How to add the nuclei to the ROI manager?

// If condition to count the ROI manager here. 

roi_count = roiManager("count"); // ROI Manager count to see the area

if (roi_count== 00) {
	print("roi < 0");
//	break; 
} 
else {
roiManager("Select", 0); // Select Nucleus in the ROI
roiManager("measure"); // Measure area

selectWindow("NR2F1"); // select the voronoi image
run("Remove Overlay"); // Remove Overlay
//roiManager("Deselect");

// Save Nuclei Measurements in Variable to save in array 
nuclei_area = getResult('Area', 0); // Area of a single cell
nuclei_mean_intensity = getResult('Mean', 0); // Mean Intensity of a single cell
nuclei_X = getResult('X', 0); // X nuclei measurement
nuclei_Y = getResult('Y', 0); // Y nuclei measurement

nr2f1_nuclei_area[i] = nuclei_area; // nr2f1 nuclei area
nr2f1_nuclei_mean_intensity[i] = nuclei_mean_intensity; // nr2f1 nuclei mean intensity
nr2d1_nuclei_x[i] = nuclei_X; // nr2f1 nuclei x coordinate
nr2f1_nuclei_y[i] = nuclei_Y; // nr2f1 nuclei  y cooreindate

wait(1000);

//if (i==10){
//break;
//	}
// Reset ROI manager 
roiManager("reset"); // Reset ROI manager 

//Array.print(nr2f1_nuclei_area); 
// Save the result in csv file

} // else
//print(nr2f1_nuclei_area[i]);
} // for loop 
//return nr2f1_nuclei_area; 


save_file_csv_loation = "D:/test.csv"; 
save_results(save_file_csv_loation,sox9_area,sox9_mean_intensity,nr2f1_nuclei_area,nr2f1_nuclei_mean_intensity);

// Function to save results to CSV file
function save_results(save_file_csv_loation, sox9_cell_area, sox9_cell_mean_intensity, nr2f1_nuclear_area, nr2f1_nuclear_mean_intensity) {
// FUNCTION TO CREATE A table to save the results
// INPUT: save_file_csv_loation, cell_area_original, cell_area_new, sox9_intensity, MenaINV_intensity
// OUTPUT: Cell Mask image title 004_voronoi_mask
// e.g. INPUT: 
// save_file_csv_loation = "C:/Users/sushukla/Desktop/temp/table2.csv";
    // cell_area_original = newArray(total_cell_counts);
    // cell_area_new = newArray(total_cell_counts);
    // sox9_intensity = newArray(total_cell_counts);
    // MenaINV_intensity = newArray(total_cell_counts);`
// e.g. OUTPUT: table2.csv


// Get the cell count array
total_cell_counts = lengthOf(sox9_cell_area); // Get the total number of cells
count_array = newArray(total_cell_counts); // Initialize the array to store the cell count
//for (i = 0; i < total_cell_counts; i++) {
for (i = 0; i < iter; i++) {

    count_array[i] = i + 1;
}

// Create table
// SOX9 --> Cell and NR2F1 --> Nuclear 
Table.create("table1"); // set a whole column
Table.setColumn("Cells_ID", count_array); // Cell ID
//Table.setColumn("SOX9_cell_area", sox9_cell_area); // Area
Table.setColumn("SOX9_cell_mean_intensity", sox9_cell_mean_intensity); // Mean Nuclear Intensity 
//Table.setColumn("NR2F1_nuclear_area", nr2f1_nuclear_area);  // Area
Table.setColumn("NR2F1_nuclear_mean_intensity", nr2f1_nuclear_mean_intensity); // Mean Nuclear Intensity 

// How to add X and Y location for cell and nuclei 

// Save Table 
saveAs("Results", save_file_csv_loation);

// Screen location
Table.setLocationAndSize(900, 500, 500, 500);
print("Results saved at: " + save_file_csv_loation);

}
run("Close All");
setBatchMode(false);