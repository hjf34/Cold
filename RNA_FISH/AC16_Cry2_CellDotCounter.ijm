///Fiji macro to count dots in cytoplasmic and nuclear area
function CytoplamNucleusCellDotCounter(x){
	
nROIs = roiManager("count");
if (nROIs > 0) {
    roiManager("Deselect"); 
	roiManager("Delete");
}

z = x + "outline.tif";

image1 = imagepath + x;
tifsav = imagepath + z;

bfwi = "open=" + image1;
run("Bio-Formats Windowless Importer", bfwi);
run("Split Channels");

chan1 = "C1-" + x;
chan2 = "C2-" + x;
chan3 = "C3-" + x;

chan3max = "MAX_" + chan3;

roisav = image1 + "roiX.zip";

///Project middle z stacks taking maximum intensity
selectWindow(chan3);
run("Z Project...", "start=9 stop=11 projection=[Max Intensity]");

////Find outline of nuclei
run("Enhance Contrast...", "saturated=3 normalize");
run("Gaussian Blur...", "sigma=2");

setAutoThreshold("Default dark");
run("Convert to Mask");
run("Fill Holes");
run("Analyze Particles...", "size=40-Infinity clear add");


nROIs = roiManager("count");

if (nROIs > 1) {
    roiManager("Combine");
    roiManager("Delete");
	roiManager("Add");
}

///required to fix bug from if only one cell
if (nROIs < 2) {
	roiManager("Select", 0);
	roiManager("Add");
	roiManager("Select", newArray(0,1));
	roiManager("Combine");
    roiManager("Delete");
	roiManager("Add");
}

/////Find outline of cells and remove holes inside cells
selectWindow(chan2);
run("Z Project...", "start=3 stop=17 projection=[Max Intensity]");
///////////////////////////////////////////////////

run("Duplicate...", " ");
setAutoThreshold("MinError");
run("Convert to Mask");
setAutoThreshold("Default dark");
run("Convert to Mask");
run("Analyze Particles...", "size=200-Infinity clear add");

nROIs = roiManager("count");

if (nROIs > 1) {
    roiManager("Combine");
    roiManager("Delete");
	roiManager("Add");
}

setAutoThreshold("Default dark");
run("Convert to Mask");
roiManager("Select", 0);
run("Make Inverse");
run("Fill Holes");
roiManager("Select", 0);
run("Make Inverse");

///////

roiManager("Show None");
roiManager("Show All");
run("Duplicate...", " ");
roiManager("Select", 0);
run("Analyze Particles...", "size=10-Infinity add");

nROIs = roiManager("count");

if (nROIs > 1) {
	roiManager("Deselect");
	roiManager("XOR");
    roiManager("Delete");
	roiManager("Add");
}

////////////////////

///////Make sure cell outline includes nuclei
selectWindow(chan3max);
roiManager("Add");
roiManager("Select", newArray(0,1));
roiManager("Combine");
roiManager("Add");
roiManager("Deselect");
roiManager("Select", 0);
roiManager("Delete");

roiManager("Save", roisav);
roiManager("Deselect"); 
roiManager("Delete"); 
run("Close All");


///////////////////////////////////////////////////

/////Creates TIF image to View Counted Dots

run("Bio-Formats Windowless Importer", bfwi);
open(roisav);

z2 = x + "CellNucCy3Dots.tif";

tifsav1 = imagepath + z2;

run("Z Project...", "start=3 stop=17 projection=[Max Intensity]");
run("Split Channels");

chan1 = "C1-MAX_" + x;
chan2 = "C2-MAX_" + x;
chan3 = "C3-MAX_" + x;

close(chan1);
close(chan3);
selectWindow(chan2);
setMinAndMax(2000, 30000);

roiManager("Select", 1);
run("Find Maxima...", "noise=6000 output=[Single Points]");
run("Analyze Particles...", "size=0-Infinity pixel add");
close();

///Outlines Nucleus and Cytoplasmic Areas
run("RGB Color");
roiManager("Show All");
setForegroundColor(0, 255, 0);
roiManager("Draw");
roiManager("Select", 0);
setForegroundColor(0, 0, 255);
roiManager("Draw");
roiManager("Select", 1);
setForegroundColor(255, 0, 0);
roiManager("Draw");

saveAs("Tiff", tifsav1);

run("Close All");

nROIs = roiManager("count");
if (nROIs > 0) {
    roiManager("Deselect"); 
	roiManager("Delete");
}

/////Image dot counting
/////Opening files

csvsav1 = image1 + "Results.csv";
bfwi = "open=" + image1;

run("Bio-Formats Windowless Importer", bfwi);
open(roisav);


//////////////////////////////////
run("Z Project...", "start=3 stop=17 projection=[Max Intensity]");
run("Split Channels");

chan1 = "C1-MAX_" + x;
chan2 = "C2-MAX_" + x;
chan3 = "C3-MAX_" + x;

close(chan1);
close(chan3);

selectWindow(chan2);
run("Clear Results");
run("Set Measurements...", "area redirect=None decimal=3");
//Nuclear area
roiManager("Select", 0);
roiManager("Measure");
//Whole cell area
roiManager("Select", 1);
roiManager("Measure");

////Cry2
///Cry2 nuclear count
roiManager("Select", 0);
run("Find Maxima...", "noise=6000 output=Count");

///Cry2 cytoplasmic count
roiManager("Select", newArray(0,1));
roiManager("XOR");
roiManager("Add");
roiManager("Select", 2);
run("Find Maxima...", "noise=6000 output=Count");

///Save Results
saveAs("Results", csvsav1);
run("Clear Results"); 

///Clear ROI
nROIs = roiManager("count");
if (nROIs > 0) {
    roiManager("Deselect"); 
	roiManager("Delete");
}

///Close All Image Files
run("Close All");

}

///Path to directory containing Cry2 Deconvolved RNA_FISH images
imagepath = "/path/to/.dv/directory/";

list1 = getFileList(imagepath);
for(i=0; i<list1.length; i++){ 
if(endsWith(list1[i], ".dv")){
CytoplamNucleusCellDotCounter(list1[i]);
}
}




