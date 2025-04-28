/* Working one!
 * Function for the manual comparison between a green/blue image and a binary image, to make side-by-side tracking more straight forward
 * Prompts the user to open two images
 * Creates an ROI around each cell of interest based upon the binary image 
 * Then goes through and zooms into the ROI's one by one, asking for user input
 * 72000
 * 20000
 * 3360
 * 
 */
 
//MAIN PROGRAM
		binaryDataArray = openBinary(); //open the binary image 
		number = binaryDataArray[0];
		binaryImageID = binaryDataArray[1];
		binaryImage = binaryDataArray[2];
		originalImageBGArray = openGBImage(); //open the GB image 
		originalImage = originalImageBGArray[0];
		originalImageID = originalImageBGArray[1];
		selectedArray = newArray(number); //array for storing comments about marked ROI's
		selectedArray = roiIterator(number, selectedArray, binaryImageID, originalImageID);
		//compiling a string with all the outputs 
		arrayString = "Index:    Comment:   \n"; //setting up the first line for output 
		for (i=0; i < number; i++) { //creating the user output to display comments about the marked microglia
			if (selectedArray[i]!=0){
				arrayString = arrayString + selectedArray[i] + "\n";
			}
		}
		//end of program message 
		Dialog.createNonBlocking("End of Program");
		Dialog.addMessage("You've hit the end of the program. Save the final edited binary image.");
		Dialog.addMessage("The following microglia ROI's were marked:");
		Dialog.addMessage(arrayString);
		Dialog.show();
		
/*HELPER FUNCTIONS
 * roiIterator(number, selectedArray, binaryImageID, originalImageID): function that iterates through ROI of interest
 * openGBImage(): function to open a GB image
 * openBinary(): function to open a binary image 
 * zoom(imageNameBinary, imageNameGB, imageIDBinary, imageIDGB, index): zooms in on two images to a highlighted ROI defined by the input index
 */
 
function roiIterator(number, selectedArray, binaryImageID, originalImageID){
/* 
 * Function to iterate through an input number of ROI's in two images
 * number: number of ROI's in the manager, integer
 * selectedArray: array for storing comments about marked microglia, array
 * binaryImageID: unique ID for the binary image 
 * originalImageID: unique ID for the GB image 
 */
	for(i=0; i < number; i++) {
		roiManager("select", i); 
		zoom(binaryImageID, originalImageID, i);
		Dialog.createNonBlocking("Manual Editor");
		Dialog.setLocation(1200,0);
		Dialog.addMessage("Make manual edits desired. Press enter to move onto next ROI");
		Dialog.addMessage("Current ROI: " + i+1);
		Dialog.addMessage("If you input f into the ROI box and press ok, the selected region in the binary image will be filled black", 8, "#FF0000");
		Dialog.addString("Enter a specific roi if desired:", "ROI"); //option to jump to a specific ROI index
		Dialog.addCheckbox("Mark Microglia", false); //option to mark for commenting
		Dialog.show();
		numberString = Dialog.getString();
		markedOrNot = Dialog.getCheckbox();
		if (markedOrNot == true) { //if you marked the ROI gives the user the ability to add a comment
			Dialog.createNonBlocking("Comments about marking");
			Dialog.addMessage("Please enter any comments about the marking, if left blank the program will default to saving only the index/name of the ROI");
			Dialog.addString("Comment: ", "Blank");
			Dialog.show();
			commentROI = Dialog.getString();
			if (commentROI == "Blank"){ //if no comment, just saves the index of the ROI
				selectedArray[i] = i+1;
			}
			else{ //if any comment saves the comment as well as the index
				selectedArray[i] = toString(i+1) + " - " + commentROI;
			}
			}
		if (numberString != "ROI"){
			if (numberString == "f"){
				selectImage(binaryImageID);
				if (selectionType != -1) {  // Check if there is a selection (ROI) active
        			setForegroundColor(0, 0, 0);  // Set the fill color to black (RGB: 0, 0, 0)
        			run("Fill");
				}
				else {
       				 print("No ROI selected");     // Inform the user if no ROI is selected
    			}
			}
			else if (isNaN(parseInt(numberString))){ //invalid input 
				Dialog.create("Erorr. Invalid ROI entered.");
				Dialog.addMessage("The program will jump to the next ROI which is: " + (toString(i+2)));
				Dialog.show();
			}
			else{ //valid input 
				i = parseFloat(numberString)-2;
			}
			}
		}
		return selectedArray;
}
function openGBImage(){
/*
 * Function to open a Green/Blue image
 * Gives the option of opening an image with RGB channels or just a regular image, only difference is whether channels tool 
 * is opened upon initiation. 
 */
		Dialog.createNonBlocking("ImageOpener");
		Dialog.addMessage("Open your original test image (blue/green)");
		Dialog.addCheckbox("Is your image a split channel RGB?", false);
		Dialog.show();
		rgbChoice = Dialog.getCheckbox();
		originalImage = File.openDialog("Open your original test image");
		open(originalImage);
		originalImageID = getImageID();
		if (rgbChoice == true){ //opening channels tool for RGB image
			run("Channels Tool...");
		}
		originalResults = newArray(2);
		originalResults[0]=originalImage;
		originalResults[1]=originalImageID;
		return originalResults;
}
function openBinary(){
/* Function that opens a user input binary image 
 * Returns an array of data from the binary image for use: 
 * [the number of ROI's in the array, the unique ID of the binary Image, the binaryImage]
 */
	//use file browser to open binary image
		print("fiji pretty please cooperate");
		Dialog.createNonBlocking("BinaryOpener");
		Dialog.addMessage("Open your binary image");
		Dialog.show();
		binaryImage = File.openDialog("Open your binary test image");
		open(binaryImage);
		binaryImageID = getImageID();
		Dialog.create("ROI Manager");
		Dialog.addCheckbox("Do you have pre-saved ROI's?", false);
		Dialog.show();
		decision = Dialog.getCheckbox();
		if (decision == true) { //if presaved ROI's
			Dialog.create("ROI opener");
			Dialog.addMessage("Choose a zip file with ROI's defined:");
			Dialog.show();
			roiOfInterest = File.openDialog("Choose your zip with ROI's defined:");
			roiManager("open", roiOfInterest);
			number = roiManager("count");		
			run("Synchronize Windows", "cursor slices frames scaling position");
		}else { //no presaved ROI's
//create an ROI around each defined area in the binary image (i.e. white sections)
			selectImage(binaryImageID);
			run("Analyze Particles...", "size=600-Infinity show=Nothing pixel include add include");
			number = roiManager("count");		
		    run("Synchronize Windows", "cursor slices frames scaling position");
	}
	resultArray = newArray(3);
	resultArray[0] = number;
	resultArray[1] = binaryImageID;
	resultArray[2] = binaryImage;
	return resultArray; 
}

function zoom(imageIDBinary, imageIDGB, index) {
/*Function to zoom to the selected ROI, takes 4 arguments
 * imageIDBinary: the unique id of the binary image, integer
 * imageIDGB: the unique id of the original image, integer
 * index: the index of the ROI that you want to zoom to, integer 
 */
// Check if there's an ROI selected
	if (roiManager("count") == 0) {
		showMessage("No ROI selected", "Please select an ROI before running this script.");
  	  	exit();
	}
// Get the selected ROI
	selectImage(imageIDBinary);
	roiManager("Select", index); // Select the first ROI (index 0)
	run("To Selection");
	selectImage(imageIDGB);
	roiManager("Select", index);
	run("To Selection");
	selectImage(imageIDGB);
	roiManager("select", 1);
}



