/*
January 2023, Qupath version compatible 0.4.1, 0.4.2
Script version 2

BSD 3-Clause License
@author Alexandre Hego
 
 contact: alexandre.hego@uliege.be
 GIGA Cell imaging facility
 University of Liege 

 **************************************************************************
 Goal:
 0) import libraries and set the variables 
 1)  Set the vectors stain (if needed)
 2) Detect Tumor annotations, merges them together 
 3) Detect if an annotation exist if not create a tissu annotation and quantify eosinophils
 3b) if a tumor annotation already exist create only a tissu annotations and continue
 4) Script to help with annotating tumor regions, separating the tumor margin from the center.
 5) remove annotation Tissu and Tumor
 6) rename annotations
 7) Detect the positive cells
 8) Save the annotations
 **************************************************************************
 
 **************************************************************************
 Tutorial
 0) use the brush tool to annotate the tumor and set it in class Tumor 
 1) Possibility to change the size of the margin for inner, outer zones
 2) Set the stainning vectors for the current batch 
 3) Detect Tumor annotations, merges them together 
 4) Detect if an annotation exist if not create a tissu annotation and quantify eosinophils
 Warning the Tumor zone need a size 2x of the expandMarginMicrons variable
 5) Select Run > Run for project
 **************************************************************************
 */


/* 0) import libraries and set the variables 
****************************************************/
import org.locationtech.jts.geom.Geometry
import qupath.lib.common.GeneralTools
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs
import static qupath.lib.gui.scripting.QPEx.*

// Set the model of random forest to detect the tissue
model_tissu = "detection_tissue"

// Set a detection for eosinophil
DetectionEosino =  runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', '{"detectionImageBrightfield":"Optical density sum","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":12.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.7,"minAreaMicrons":15.0,"maxAreaMicrons":350.0,"threshold":0.1,"maxBackground":2.0,"watershedPostProcess":true,"excludeDAB":false,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true,"thresholdCompartment":"Cell: DAB OD mean","thresholdPositive1":0.2,"thresholdPositive2":0.4,"thresholdPositive3":0.6000000000000001,"singleThreshold":true}');


// How much to expand each region
double expandMarginMicrons = 1200

// Define the colors
def colorInnerMargin = getColorRGB(0, 0, 200)
def colorOuterMargin = getColorRGB(0, 200, 0)
def colorCentral = getColorRGB(200, 0, 0)

/* 1)  Set the vectors stain (if needed)
****************************************************/
//setColorDeconvolutionStains('{"Name" : "H-DAB estimated", "Stain 1" : "Hematoxylin", "Values 1" : "0.66034 0.69836 0.27614", "Stain 2" : "DAB", "Values 2" : "0.50613 0.58015 0.63817", "Background" : " 187 182 174"}');

/* 2) Detect Tumor annotations, merges them together 
****************************************************/
selectObjectsByClassification("Tumor");
mergeSelectedAnnotations();
resetSelection();

/* 3) Detect if an annotation exist if not create a tissu annotation and quantify eosinophils
*********************************************************************************************/
count = 0
for (annotation in getAnnotationObjects()) {
    count = count +1
}

if (count == 0) {
    createAnnotationsFromPixelClassifier(model_tissu, 5000000.0, 0.0)
    Thread.sleep(10)
    selectObjectsByClassification("Tissu")
    DetectionEosino

    path = buildFilePath(PROJECT_BASE_DIR, 'Measurements')
    name = getProjectEntry().getImageName() + '.tsv'
    //make sure the directory exists
    mkdirs(path)
    // Save the results
    path = buildFilePath(path, name)
    selectObjectsByClassification("Tissu")
    saveAnnotationMeasurements(path)
} 

else {
        createAnnotationsFromPixelClassifier(model_tissu, 5000000.0, 0.0)

/* 4) Script to help with annotating tumor regions, separating the tumor margin from the center.
***************************************************************************************************/
// Choose whether to lock the annotations or not (it's generally a good idea to avoid accidentally moving them)
def lockAnnotations = true

// Extract the main info we need
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()

// We need the pixel size
def cal = server.getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) {
  print 'We need the pixel size information here!'
  return
}
if (!GeneralTools.almostTheSame(cal.getPixelWidthMicrons(), cal.getPixelHeightMicrons(), 0.0001)) {
  print 'Warning! The pixel width & height are different; the average of both will be used'
}

// Get annotation & detections
hierarchy = getCurrentHierarchy()
annotation = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("Tumor")}
hierarchy.getSelectionModel().setSelectedObject(annotation)

def annotations = getAnnotationObjects()
def selected = getSelectedObject()
if (selected == null || !selected.isAnnotation()) {
  print 'Please select an annotation object!'
  return
}

// We need one selected annotation as a starting point; if we have other annotations, they will constrain the output
annotations.remove(selected)

// Extract the ROI & plane
def roiOriginal = selected.getROI()
def plane = roiOriginal.getImagePlane()

// If we have at most one other annotation, it represents the tissue
Geometry areaTissue
PathObject tissueAnnotation
if (annotations.isEmpty()) {
  areaTissue = ROIs.createRectangleROI(0, 0, server.getWidth(), server.getHeight(), plane).getGeometry()
} else if (annotations.size() == 1) {
  tissueAnnotation = annotations.get(0)
  areaTissue = tissueAnnotation.getROI().getGeometry()
} else {
  print 'Sorry, this script only support one selected annotation for the tumor region, and at most one other annotation to constrain the expansion'
  return
}

// Calculate how much to expand the middle layer
double expandPixels = expandMarginMicrons / cal.getAveragedPixelSizeMicrons()
def areaTumor = roiOriginal.getGeometry()

// Get the outer margin area
def geomInner = areaTumor.buffer(expandPixels)
geomInner = geomInner.difference(areaTumor)
geomInner = geomInner.intersection(areaTissue)
def roiInner = GeometryTools.geometryToROI(geomInner, plane)
def annotationInner = PathObjects.createAnnotationObject(roiInner, getPathClass("Inner"))
annotationInner.setName("Inner margin")
annotationInner.setColorRGB(colorInnerMargin)
addObject(annotationInner)

// Calculate how much to expand the outer layer
double expandPixels2 = expandMarginMicrons*2 / cal.getAveragedPixelSizeMicrons()

// Get the outer margin area
def geomOuter = areaTumor.buffer(expandPixels2)
geomOuter = geomOuter.difference(geomInner)
geomOuter = geomOuter.difference(areaTumor)
geomOuter = geomOuter.intersection(areaTissue)
def roiOuter = GeometryTools.geometryToROI(geomOuter, plane)
// def annotationOuter = PathObjects.createAnnotationObject(roiOuter)
def annotationOuter = PathObjects.createAnnotationObject(roiOuter, getPathClass("Outer"))
annotationOuter.setName("Outer margin")
annotationOuter.setColorRGB(colorOuterMargin)
addObject(annotationOuter)

// Get the central area
geomCentral = areaTumor.intersection(areaTissue)
def roiCentral = GeometryTools.geometryToROI(geomCentral, plane)
// change def annotationCentral = PathObjects.createAnnotationObject(roiCentral)
def annotationCentral = PathObjects.createAnnotationObject(roiCentral, getPathClass("Center"))
annotationCentral.setName("Center")
annotationCentral.setColorRGB(colorCentral)
addObject(annotationCentral)

/* 5) remove annotation Tissu and Tumor
***************************************************************************************************/

selectObjectsByClassification("Tissu")
clearSelectedObjects(false)

selectObjectsByClassification("Tumor")
clearSelectedObjects(false)



/* 6) rename annotations
***************************************************************************************************/
for (annotation in getAnnotationObjects()) {
    if (annotation.getPathClass() == getPathClass("Outer")) {
    annotation.setPathClass(getPathClass("Outer_external")) 
    }
    if (annotation.getPathClass() == getPathClass("Inner")) {
    annotation.setPathClass(getPathClass("Outer_internal")) 
    }
    if (annotation.getPathClass() == getPathClass("Center")) {
    annotation.setPathClass(getPathClass("Tumor")) 
    } 
}


/* 7) Detect the positive cells
***************************************************************************************************/
selectObjectsByClassification("Outer_external")
runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', '{"detectionImageBrightfield":"Optical density sum","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":12.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.7,"minAreaMicrons":15.0,"maxAreaMicrons":350.0,"threshold":0.1,"maxBackground":2.0,"watershedPostProcess":true,"excludeDAB":false,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true,"thresholdCompartment":"Cell: DAB OD mean","thresholdPositive1":0.2,"thresholdPositive2":0.4,"thresholdPositive3":0.6000000000000001,"singleThreshold":true}');

selectObjectsByClassification("Outer_internal")
runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', '{"detectionImageBrightfield":"Optical density sum","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":12.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.7,"minAreaMicrons":15.0,"maxAreaMicrons":350.0,"threshold":0.1,"maxBackground":2.0,"watershedPostProcess":true,"excludeDAB":false,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true,"thresholdCompartment":"Cell: DAB OD mean","thresholdPositive1":0.2,"thresholdPositive2":0.4,"thresholdPositive3":0.6000000000000001,"singleThreshold":true}');

selectObjectsByClassification("Tumor")
runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', '{"detectionImageBrightfield":"Optical density sum","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":12.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.7,"minAreaMicrons":15.0,"maxAreaMicrons":350.0,"threshold":0.1,"maxBackground":2.0,"watershedPostProcess":true,"excludeDAB":false,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true,"thresholdCompartment":"Cell: DAB OD mean","thresholdPositive1":0.2,"thresholdPositive2":0.4,"thresholdPositive3":0.6000000000000001,"singleThreshold":true}');



/* 8) Save the annotations
***************************************************************************************************/
path = buildFilePath(PROJECT_BASE_DIR, 'Measurements')

name = getProjectEntry().getImageName() + '.tsv'


//make sure the directory exists
mkdirs(path)

// Save the results
path = buildFilePath(path, name)
selectObjectsByClassification("Tumor")
saveAnnotationMeasurements(path)

}