// This script was written by Friederike Kessel,
// In the Workgroup of Prof. Hugo at the University Hospital Carl Gustav Carus
// Of the medical faculty of the University Dresden

// Working Directory
argument=getArgument();
path=argument+"\\"

// Setup
run("Bio-Formats Macro Extensions");
close("*");
roiManager("Reset");
roiManager("Show All");
run("Set Measurements...", "standard median area_fraction redirect=None decimal=6");
run("Clear Results");

// Files
File.makeDirectory(path+"Results");
files=getFileList(path);

// PART 1: SETTINGS--------------------------------------------------------------------
run("Clear Results");

for(a=0; a<files.length;a++){
	print(files[a]);
	if(endsWith(files[a], ".lif")){		
		Ext.setId(path+files[a]);
		Ext.getSeriesCount(seriesCount);		
		for(b=0; b<seriesCount; b++){
			Ext.setSeries(b);
			Ext.getSeriesName(seriesName);			
			if(matches(seriesName, ".*LY.*")){
				print(seriesName);
				name1=substring(files[a], 0, lengthOf(files[a])-4);				
				setResult("File", nResults, name1);
				setResult("Series_LY", nResults-1, seriesName);
				setResult("Series#_LY", nResults-1, b);							
				run("Bio-Formats Importer", "open=["+path+files[a]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+b+1);				
				rename("Current");
					
				settingsIntensity();

				// To find the position of the time series in the z-stack later
				selectWindow("Current");
				getDimensions(width, height, channels, slices, frames);
				setSlice(slices);
				run("Duplicate...", "slice");
				rename("Compare");
				close("Current");							
				roiManager("Save", path+"Results\\"+name1+"_"+seriesName+".zip");

				//Manually selection corresponding z-stack
				Dialog.create("Select corresponding z-stack");
				Dialog.addMessage("For sample "+seriesName);
				
				serieslist=newArray(1);
				for(c=0; c<seriesCount; c++){
					Ext.setSeries(c);
					Ext.getSeriesName(seriesName);
					if(matches(seriesName, ".*PT.*")){
						serieslist=Array.concat(serieslist, seriesName);
					}
				}
				
				Dialog.addChoice("Series", serieslist);				
				Dialog.show();
				zseries=Dialog.getChoice();

				// settings for volume of the sample
				setBatchMode(true);
				for(c=0; c<seriesCount; c++){
					Ext.setSeries(c);
					Ext.getSeriesName(seriesName);					
					if(seriesName==zseries){
						if(matches(seriesName, ".*PT.*")){						
							setResult("Series_PT", nResults-1, zseries);
							setResult("Series#_PT", nResults-1, c);
							updateResults();													
							setResult("low_limit", nResults-1, 0);
							setResult("high_limit", nResults-1, 100);
							updateResults();							
							roiManager("Save", path+"Results\\Volume_"+name1+"_"+seriesName+".zip");
						}
					}
				}
				setBatchMode(false);
				roiManager("reset");	
				close("*");
			}			
		}
	}
}


function settingsIntensity(){
	selectImage("Current");
	run("Duplicate...", "duplicate channels=2-2");
	close("Current");
	rename("Current");

	// Preprocess stack-----------------------------------------------------
	setLocation(0,0,800,800);	
	run("Median...", "radius=1 stack");
	run("Enhance Contrast...", "saturated=0");

	getDimensions(width, height, channels, slices, frames);
	setSlice(frames/4*3);
			
	// line along center of tubule
	setTool("polyline");
	run("Line Width...", "line=10");
	waitForUser("Draw a line along the center of the tubule");
	while(selectionType()==-1){
		waitForUser("Draw a line along the center of the tubule");
	}
	
	run("Fit Spline");
	roiManager("Add");
	roiManager("Select", 0);
	roiManager("remove slice info");	
}



// save settings
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
saveAs("Results", path+"Results\\Settings_"+year+"_"+month+"_"+dayOfMonth+"_"+hour+"-"+minute+"-"+second+".txt");
		
files=newArray(nResults);
Series_LY=newArray(nResults);
SeriesNo_LY=newArray(nResults);
Series_PT=newArray(nResults);
SeriesNo_PT=newArray(nResults);
low_limit=newArray(nResults);
high_limit=newArray(nResults);

for(a=0; a<nResults; a++){
	files[a]=getResultString("File", a);
	Series_LY[a]=getResultString("Series_LY", a);
	SeriesNo_LY[a]=getResult("Series#_LY", a);
	Series_PT[a]=getResultString("Series_PT", a);
	SeriesNo_PT[a]=getResult("Series#_PT", a);
	low_limit[a]=getResult("low_limit", a);
	high_limit[a]=getResult("high_limit", a);
}

run("Clear Results");

// PART 2 MEASUREMENT
setBatchMode(true);
for(a=0; a<files.length;a++){
	print(files[a]);	
	Ext.setId(path+files[a]+".lif");
	Ext.getSeriesCount(seriesCount);		
	run("Bio-Formats Importer", "open=["+path+files[a]+".lif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+SeriesNo_LY[a]+1);					
	rename("Current");					
	roiManager("Open", path+"Results\\"+files[a]+"_"+Series_LY[a]+".zip");	
	
	measureIntensity();
	
	close("Current");
	selectWindow("Nuclei");
	getDimensions(width, height, channels, slices, frames);
	setSlice(slices);
	run("Duplicate...", "slice");
	rename("Compare");
	close("Nuclei");
				
	saveAs("Results", path+"Results\\"+files[a]+"_"+Series_LY[a]+".txt");
	run("Clear Results");
	roiManager("Reset");


	//measure volume of sample
	if(SeriesNo_PT[a]!=0){
		Ext.setSeries(SeriesNo_PT[a]);

		run("Bio-Formats Importer", "open=["+path+files[a]+".lif] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+SeriesNo_PT[a]+1);
		rename("Volume");
		getVoxelSize(width, height, depth, unit);
		vwidth=width;
		vdepth=depth;
		roiManager("Open", path+"Results\\Volume_"+files[a]+"_"+Series_PT[a]+".zip");	
		
		preprocess_stack();

		compare_stack_position()


		
		roiManager("Save", path+"Results\\Volume_"+files[a]+"_"+Series_PT[a]+".zip");
		filepath=path+"Results\\Volume_"+files[a]+"_"+Series_LY[a];
								
		measureVolume(filepath);
		updateResults();
		saveAs("Results", path+"Results\\Volume_"+files[a]+"_"+Series_LY[a]+".txt");
		roiManager("reset");	
		close("*");		
	}
}


function measureIntensity(){
	selectImage("Current");
	run("Duplicate...", "duplicate channels=4-4");
	rename("Nuclei");
	selectImage("Current");
	run("Duplicate...", "duplicate channels=2-2");
	close("Current");
	rename("Current");
	
	//Preprocess stack-----------------------------------------------------
	setLocation(500,20,1000,1000);	
	run("Median...", "radius=1 stack");
	run("Enhance Contrast...", "saturated=0");
	getDimensions(width, height, channels, slices, frames);
	setSlice(frames/4*3);
			
	//line along center of tubule
	roiManager("Select", 0);
	run("Plots...", "width=530 height=300 font=12 draw draw_ticks auto-close minimum=0 maximum=0 interpolate");
	getDimensions(width, height, channels, slices, frames);
	run("Clear Results");
	for(c=1; c<frames+1; c++){	
		selectWindow("Current");
		setSlice(c);
		run("Plot Profile");
		Plot.getValues(x, y);
					 
		 for (i=0; i<x.length; i++){
		     setResult("length", nResults, x[i]);
		     setResult("intensity", nResults-1, y[i]);
		     setResult("frame", nResults-1, c);
		 }
 		close();
			
	}
}

function preprocess_stack(){
	//Preprocess stack-----------------------------------------------------
	run("Split Channels");
	close("C1-Volume");
	getDimensions(width, height, channels, slices, frames);
	setSlice(slices/1.5);
	imglist=getList("image.titles");
	getVoxelSize(width, height, depth, unit);
	vwidth=width;
	vdepth=depth;
			
	for(c=0;c<imglist.length; c++){
		selectWindow(imglist[c]);
		if(imglist[c]!="Compare"){
			setSlice(slices/1.5);
			run("Median 3D...", "x="+3*vwidth+" y="+3*vwidth+" z="+3*vdepth);
			run("Enhance Contrast", "saturated=1");
			run("Apply LUT", "stack");
		}
	}
}

function compare_stack_position(){
	selectImage("Compare");
	getVoxelSize(width, height, depth, unit);
	vwidth=width;
	run("Median...", "radius="+3*vwidth);
	roiManager("Select", 0);
	run("Line to Area");
	run("Enlarge...", "enlarge=15 pixel");
	run("Enhance Contrast", "saturated=1");
	run("Select None");
	run("Apply LUT");	
	imageCalculator("Difference create stack", "C4-Volume","Compare");
	close("Compare");
	rename("Compare");
	run("Invert", "stack");
		
	getDimensions(width, height, channels, slices, frames);
	xmean=newArray(slices);
	xsd=newArray(slices);
	xcompare=newArray(slices);
	targetslice=1;
	finalcompare=0;
	for(d=1; d<=slices; d++){
		setSlice(d);
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		xmean[d-1]=mean;
						
		if(xmean[d-1]>finalcompare){
			finalcompare=xmean[d-1];
			targetslice=d;
		}
	}
	print(targetslice);
	selectWindow("C2-Volume");
	getDimensions(width, height, channels, slices, frames);
	setSlice(targetslice);
	roiManager("Select", 0);
	setSlice(targetslice);
	roiManager("Update");
}
		
function measureVolume(filepath){	
	selectImage("C4-Volume");
	getVoxelSize(width, height, depth, unit);
	vwidth=width;
	vdepth=depth;
	run("Variance...", "radius=2 stack");
	run("Maximum 3D...", "x="+4*vwidth+" y="+4*vwidth+" z="+4*vdepth);
	run("Median 3D...", "x="+4*vwidth+" y="+4*vwidth+" z="+4*vdepth);
	run("Convert to Mask", "method=Default background=Dark calculate black");
	
	selectImage("C2-Volume");
	run("Subtract Background...", "rolling=50 stack");
	imageCalculator("Subtract create stack", "C2-Volume","C3-Volume");
	imageCalculator("Subtract create stack", "Result of C2-Volume","C4-Volume");
	rename("Volume");
	
	close("C2-Volume");
	close("C3-Volume");
	close("C4-Volume");
	close("Result of C2-Volume");	
	setLocation(500,20,1000,1000);

	//Mask ROI-----------------------------------------------------
	getDimensions(width, height, channels, slices, frames);
	setSlice(slices/1.5);				
	roiManager("Select", 0);
	run("Line to Area");
	run("Enlarge...", "enlarge=15 pixel");
	run("Clear Outside", "stack");
	run("Enhance Contrast...", "saturated=0");
	run("Apply LUT", "stack");
	run("Select None");
	selectImage("Volume");	
		
	//Measure intensity					
	roiManager("Select", 0);
	Roi.getPosition(channel, slice, frame);
	sliceno=slice;
	run("Line to Area");
	roiManager("Add");
	setThreshold(1, 255);
	run("Create Selection");
	if(selectionType()==-1){
		setBatchMode("exit and display");
		exit("Measure intensity	in volume measurement: no threshold");
	}
	roiManager("Add");

	roiManager("Select", newArray(1,2));
	roiManager("AND");
	getStatistics(area, mean, min, max, std, histogram);
	run("Select None");
	roiManager("Select", newArray(1,2));
	roiManager("Delete");
	run("Select None");				
					
	//Reduce z-stack
	selectImage("Volume");
	run("Select None");
	rename("ZReduced");

	//Clear Outside selection
	getDimensions(width, height, channels, slices, frames);
	setBackgroundColor(0,0,0);
	for(c=1; c<=slices; c++){
		setSlice(c);
		setThreshold(2*mean, 255);
		run("Create Selection");
		
		if(selectionType()!=-1){
			run("Clear", "slice");
		}
	}

	run("Select None");									
	run("Scale...", "x=0.5 y=0.5 z=1.0 interpolation=Bilinear average process create");
	rename("Scaled");
	close("ZReduced");

	selectImage("Scaled");
	roiManager("Select", 0);		
	run("Scale... ", "x=0.5 y=0.5");
	getSelectionCoordinates(xpoints, ypoints);
	run("Select None");
	run("Duplicate...", "duplicate");
	setSlice(sliceno);
	makeRectangle(xpoints[xpoints.length/2], ypoints[ypoints.length/2], 1,1);
	run("Fill", "slice");
	rename("peaks");
	
	//Watershed

	run("3D Watershed", "seeds_threshold=1 image_threshold="+mean-std+" image=Scaled seeds=peaks radius=20");
	run("8-bit");
	run("Select None");
	run("Scale...", "x=2 y=2 z=1.0 interpolation=Bilinear average process create");
	close("watershed");
	rename("watershed");
	close("peaks");
	
	//Convert stack to mask					
	getDimensions(width, height, channels, slices, frames);
	setForegroundColor(255,255,255);
	for(c=1; c<=slices; c++){
		setSlice(c);
		setThreshold(1, 65535);
		run("Create Selection");
		if(selectionType()!=-1){
			run("Fill", "slice");
		}
	}
	run("Remove Outliers...", "radius=20 threshold=1 which=Bright stack");
	run("Remove Outliers...", "radius=30 threshold=100 which=Dark stack");	
	run("Median 3D...", "x="+4*vwidth+" y="+4*vwidth+" z="+4*vdepth);			
	run("Select None");

	//Find biggest volume,	
	run("3D Manager Options", "volume distance_between_centers=5 distance_max_contact=0.9 drawing=Contour");
	run("Clear Results");
	run("3D Simple Segmentation", "low_threshold=1 min_size=0 max_size=-1");
	close("Bin");
	close("watershed");
	run("3D Geometrical Measure");
	max=0;
	maxobj=0;
	for(b=0; b<nResults; b++){
		V=getResult("Volume(pix)", b);
		if(V>max){
			max=V;
			maxobj=b;
			value=getResult("Value", b);
		}
	}
	run("Colors...", "foreground=white background=black selection=magenta");
	if(max>0){
		//erase smaller blobs
		getDimensions(width, height, channels, slices, frames);
		for(c=1; c<=slices; c++){
			setSlice(c);
			setThreshold(value, value);
			run("Create Selection");
			if(selectionType()!=-1){
				setBackgroundColor(0, 0, 0);
				run("Clear Outside", "slice");
				setForegroundColor(255,255,255);
				run("Fill", "slice");
				run("Select None");

			}else{
				run("Select All");
				run("Clear", "slice");
			}
		}
	}

	run("Maximum 3D...", "x="+3*vwidth+" y="+3*vwidth+" z="+3*vdepth);
	run("Remove Outliers...", "radius="+10*vwidth+" threshold=1 which=Bright stack");
	run("8-bit");

	run("Properties...", "unit=micron pixel_width="+vwidth+" pixel_height="+vwidth+" voxel_depth="+vdepth);
	getDimensions(width, height, channels, slices, frames);
		
	saveAs("TIFF", filepath);
	rename("watershed");
	run("Clear Results");
			
	//Get sample volume, slice by slice (Measurement of tubular volume up to each position)
	roiManager("Select", 0);
	run("Fit Spline");
	roiManager("Update");
					
	getSelectionCoordinates(x, y);
	volume=newArray(x.length-1);
	
	for(c=0; c<x.length-1; c++){
		run("Select None");
		run("Duplicate...", "duplicate");
		makeLine(x[c], y[c], x[c+1], y[c+1]);
		run("Rotate...", "  angle=90");
		getLine(x1, y1, x2, y2, lineWidth);
		is_length=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
						
		delta_x=(x1-x2)*50*vwidth/is_length;
		delta_y=(y1-y2)*50*vwidth/is_length;
				
		x1_x=x1-((delta_x/2)-(x1-x2));
		x2_x=x2+((delta_x/2)-(x1-x2));
		y1_x=y1-((delta_y/2)-(y1-y2));
		y2_x=y2+((delta_y/2)-(y1-y2));
		makeLine(x1_x, y1_x, x2_x, y2_x);
		run("Line to Area");
		run("Enlarge...", "enlarge=1 pixel");
		roiManager("Add");
		if(c>0){
			roiManager("Select", newArray(roiManager("count")-1, roiManager("count")-2));
			roiManager("Combine");
						
			roiManager("Add");
			roiManager("Select", roiManager("count")-2);
			roiManager("Delete");
		}
		roiManager("Select", roiManager("count")-1);
		run("Clear Outside", "stack");
		run("Select None");
					
		run("3D Geometrical Measure");
		volumex=0;
		for(d=0; d<nResults; d++){
			volumex=volumex+getResult("Volume(unit)", d);
		}
		volume[c]=volumex;
		open_img=getList("image.titles");
		for(d=0; d<open_img.length; d++){
			if(open_img[d]!="watershed"){
				close(open_img[d]);
			}
		}
		run("Clear Results");
	}

	for(c=0; c<x.length-1; c++){
		setResult("Position", nResults, c+1);
		setResult("Volume", nResults-1, volume[c]);
	}
}
setBatchMode(false);