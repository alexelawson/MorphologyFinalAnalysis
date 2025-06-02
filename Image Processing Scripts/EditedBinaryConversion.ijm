/*
 * ImageJ Macro to take an RGB image (or folder with RGB images) and convert them to binary. 
 * Currently with a pre-set threshold: 50, 255
 * Created by Alex Lawson 
 */
 
Dialog.create("Choice");
Dialog.addMessage("Open the folder with images you want to convert to binary.");
Dialog.addCheckbox("Check here if only interested in processing 1 image", false);
Dialog.show();
userInputChoice = Dialog.getCheckbox();

if (userInputChoice==true){ //if user only wants to open one image
	singleBinaryConversion();
}
else{ //if user wants to open multiple images 
	
	//opens the file browser for user to choose a directory with images 
	setOption("JFileChooser",true);
	image_dir=getDirectory("Choose parent folder containing images");
	image_input=getFileList(image_dir); //creating an array of the images that you want to process 
	count=image_input.length;
		
	//opens the file browser for user to choose a directory to save images to 
	setOption("JFileChooser",true);
	binary_output=getDirectory("Choose output folder to write binary images to");
	
	//iterate through each image, convert to binary, and save in the selected directory 
	for (i=(0); i<(count); i++){ 
		multipleBinaryConversion(image_dir, binary_output, image_input[i]);
	}
}
Dialog.create("End of Program");
Dialog.addMessage("Program completed. \n If you processed a folder with images they should be saved in your selected directory. \n If you only processed 1 image, the program DID NOT save the image and you can do so now.");
Dialog.show();


/*
 * Function to convert an image into binary 
 * Args: img_name = string that's the name of the image you want to process 
 */
function imageConversion(img_name){
	img_namenoext = replace(img_name , ".tif" , "" );
	run("Make Composite");
	run("Split Channels");
	selectWindow("C2-"+img_name);
	run("Duplicate...", " ");
	selectWindow("C3-"+img_name);
	run("Close");
	selectWindow("C1-"+img_name);
	run("Close");
	selectWindow("C2-"+img_name);
	run("Close");
	selectWindow("C2-"+img_namenoext+"-1.tif");
	run("8-bit");
	//convert to grayscale 
	run("Grays");
	run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.35");
	//threshold set to retain the most possible microglia, without adding detail where there is none
	//autothresholds could also be used, but as we are manually editing the data, we wanted to include
	//rather than exclude sections
	setThreshold(50, 255);	//modify as eeded and if it works for your images
	run("Despeckle");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Dark dark");
	run("Close-"); // this function connects two dark pixels if they are separated by up to 2 pixels (helpful for skeletonization)
	//replaces a bright or dark outlier pixel by the median of the pixels in the surrounding area	//area is set as a radius of 2, threshold set to define an outlier as anything >50% different
	run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
	//Removes any cells below 600pix area using the white cell mask we created earlier (you can modify this as needed)
	run("Analyze Particles...", "size=600-Infinity pixel show=[Masks]"); //this step is unnecessary as if you run the ciernia lab's protocol you will specify this in the single cell conversion, however we included with our preset for testing purposes
	run("Invert LUTs");
	selectWindow("C2-"+img_namenoext+"-1.tif");
	run("Close");
}

//Function to convert multiple images to binary 
// ARGs: input = a directory with images you want to convert, 
//       output = a directory path where you want to save images
//       filename = the specific image that you want to process  
function multipleBinaryConversion(input, output, filename){
	print(input + filename);
    open(input + filename);
    img_name=getTitle();
	dirCropOutput=output;
	imageConversion(img_name);
	notiff = replace(img_name, ".tif", "");
	saveAs("Tiff", dirCropOutput+ notiff);			
	selectWindow(img_name);
	run("Close");
    }
    
//Function to convert a single image to binary     
function singleBinaryConversion(){
	Dialog.create("ImageOpener");
	Dialog.addMessage("Open your single image.");
	Dialog.show();
	originalImage = File.openDialog("Open your original test image");
	open(originalImage);
	originalImageID = getImageID();
	img_name = getTitle();
	imageConversion(img_name);
}
