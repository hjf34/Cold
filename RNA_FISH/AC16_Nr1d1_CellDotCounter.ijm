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

chan2max = "MAX_" + chan2;

roisav = image1 + "roiX.zip";

///Project middle z stacks taking maximum intensity
selectWindow(chan2);
run("Z Project...", "start=14 stop=17 projection=[Max Intensity]");

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
selectWindow(chan1);
run("Z Project...", "start=5 stop=25 projection=[Max Intensity]");
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
selectWindow(chan2max);
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

/////Remove very bright dots/areas from Nucleus (Non specific probe binding)

run("Bio-Formats Windowless Importer", bfwi);

run("Z Project...", "start=5 stop=25 projection=[Max Intensity]");
run("Split Channels");

chan1 = "C1-MAX_" + x;
chan2 = "C2-MAX_" + x;

close(chan2);
run("Duplicate...", " ");

setThreshold(4000, 65535);
run("Convert to Mask");
run("Analyze Particles...", "size=8-Infinity pixel clear add");

nROIs = roiManager("count");

if(nROIs>0){
	
	if (nROIs > 1) {
		roiManager("Combine");
    	roiManager("Delete");
		roiManager("Add");
    }

	open(roisav);

	roiManager("Select", newArray(0,1));
	roiManager("AND");
	type1 = selectionType();
	if (type1==9) {
		roiManager("Add");
		roiManager("Select", newArray(0,2));
		roiManager("AND");
		roiManager("Add");
		roiManager("Select", newArray(1,3));
		roiManager("XOR");
		roiManager("Add");
		roiManager("Select", newArray(2,4));
		roiManager("XOR");
		roiManager("Add");
		roiManager("Select", newArray(0,1,2,3,4));
		roiManager("Delete");
	}

///if no bright dots in nuclei
	if (type1==-1) {
		roiManager("Select", newArray(0,2));
		roiManager("AND");
		type2 = selectionType();
		if (type2==9){
			roiManager("Add");
			roiManager("Select", newArray(2,3));
			roiManager("XOR");
			roiManager("Add");
			roiManager("Select", newArray(0,2,3));
			roiManager("Delete");
		}
		if (type2==-1){
			roiManager("Deselect");
			roiManager("Select", 0);
			roiManager("Delete");
		}
}
/////
roiManager("Save", roisav);
}
	

nROIs = roiManager("count");
if (nROIs > 0) {
    roiManager("Deselect"); 
	roiManager("Delete");
}

run("Close All");

////////////////////////////////////////////////
////////////////////////////////////////////////

///////////////////////////////////////////////////
///////////////////////////////////////////////////

/////Creates TIF image to View Counted Dots

run("Bio-Formats Windowless Importer", bfwi);
open(roisav);

z2 = x + "CellNucCy3Dots.tif";

tifsav1 = imagepath + z2;

run("Z Project...", "start=5 stop=25 projection=[Max Intensity]");
run("Split Channels");

chan1 = "C1-MAX_" + x;
chan2 = "C2-MAX_" + x;

close(chan2);
selectWindow(chan1);
setMinAndMax(1000, 6000);

roiManager("Select", 1);
run("Find Maxima...", "noise=700 output=[Single Points]");
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
run("Z Project...", "start=5 stop=25 projection=[Max Intensity]");
run("Split Channels");

chan1 = "C1-MAX_" + x;
chan2 = "C2-MAX_" + x;

selectWindow(chan1);
run("Clear Results");
run("Set Measurements...", "area redirect=None decimal=3");
//Nuclear area
roiManager("Select", 0);
roiManager("Measure");
//Whole cell area
roiManager("Select", 1);
roiManager("Measure");

///NR1D1 nuclear count
roiManager("Select", 0);
run("Find Maxima...", "noise=700 output=Count");

///NR1D1 cytoplasmic count
roiManager("Select", newArray(0,1));
roiManager("XOR");
roiManager("Add");
roiManager("Select", 2);
run("Find Maxima...", "noise=700 output=Count");

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

///Path to directory containing Nr1d1 Deconvolved RNA_FISH images (AC16 cells)
imagepath = "/path/to/.dv/directory/";

list1 = getFileList(imagepath);
for(i=0; i<list1.length; i++){ 
if(endsWith(list1[i], ".dv")){
CytoplamNucleusCellDotCounter(list1[i]);
}
}




