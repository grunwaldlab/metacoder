/**
 *
 * @author Yoshiki Vazquez Baeza
 * @copyright Copyright 2013, The Emperor Project
 * @credits Yoshiki Vazquez Baeza
 * @license BSD
 * @version 0.9.5
 * @maintainer Yoshiki Vazquez Baeza
 * @email yoshiki89@gmail.com
 * @status Release
 *
 */

/**
 *
 * @name TrajectoryOfSamples
 *
 * @class Represents an ordered set of samples where each sample is indexed by
 * the sample identifier and the order is as given by their position along the
 * gradient.
 *
 * @property {Array} [sampleNames=Array()] array of sample identifiers.
 * @property {Array} [gradientPoints=Array()] array of values where each of the
 * samples exist in.
 * @property {float} [minimumDelta=float] minimum differential between samples
 * in the trajectory; the value is computed using the gradientPoints array.
 * @property {float} [suppliedN=float] minimum number of frames a distance will
 * have in the gradient.
 * @property {Array} [interpolatedCoordinates=Array()] array of objects with
 * the corresponding interpolated x, y and z values. The interpolation
 * operation takes place between subsequent samples.
 *
 */


/**
 *
 * @name TrajectoryOfSamples
 *
 * @class Represents an ordered set of samples and their position in PCoA
 * space.
 *
 * @param {sampleNames} an Array of strings where each string is a sample
 * identifier.
 * @param {metadataCategoryName} a string indicating the name of the category
 * in the mapping file used to generate this trajectory.
 * @param {gradientPoints} an Array of floating point values where each value
 * corresponds to the position of the samples in the gradient.
 * @param {coordinates} an Array of objects with x, y and z properties where
 * each corresponds to the position of a sample in PCoA space.
 * @param {minimumDelta} minimum differential between the ordered
 * gradientPoints this value must be non-zero. Note that this value should be
 * computed taking into account all the other trajectories that will be
 * animated together, usually by an AnimationDirector object.
 * @param {suppliedN} a parameter that will determine how many points should
 * should be found in the the trajectory.
 * @param {maxN} maximum number of samples allowed per interpolation interval.
 *
 **/
function TrajectoryOfSamples(sampleNames, metadataCategoryName, gradientPoints,
                             coordinates, minimumDelta, suppliedN, maxN){
  this.sampleNames = sampleNames;
  this.metadataCategoryName = metadataCategoryName;

  // array of the values that samples have through the gradient
  this.gradientPoints = gradientPoints;

  // the first three axes of the data points
  this.coordinates = coordinates;

  // minimum distance in the gradient
  this.minimumDelta = minimumDelta;

  // this value determines how fast the animation will run for now let's use
  // 5 and stick to it as a good default value; 60 was way too slow
  this.suppliedN = suppliedN !== undefined ? suppliedN : 5;
  this.maxN = maxN !== undefined ? maxN : 10;

  if (this.coordinates.length != this.gradientPoints.length) {
    throw new Error("The number of coordinate points and gradient points is"+
                    "different, make sure these values are consistent.");
  }

  // initialize as an empty array but fill it up upon request
  this.interpolatedCoordinates = null;
  this._generateInterpolatedCoordinates();

  return this;
}

/**
 *
 * Helper method to iterate over all the coordinates and generate interpolated
 * arrays.
 *
 */
TrajectoryOfSamples.prototype._generateInterpolatedCoordinates = function(){
  var pointsPerStep = 0, delta = 0;
  var interpolatedCoordinatesBuffer = new Array(),
      intervalBuffer = new Array();
  var currInterpolation;

  // iterate over the gradient points to compute the interpolated distances
  for (var index = 0; index < this.gradientPoints.length-1; index++){

    // calculate the absolute difference of the current pair of points
    delta = Math.abs(Math.abs(this.gradientPoints[index])-Math.abs(
        this.gradientPoints[index+1]));

    pointsPerStep = this.calculateNumberOfPointsForDelta(delta);
    if (pointsPerStep > this.maxN){
      pointsPerStep = this.maxN;
    }

    currInterpolation = linearInterpolation(this.coordinates[index]['x'],
                                            this.coordinates[index]['y'],
                                            this.coordinates[index]['z'],
                                            this.coordinates[index+1]['x'],
                                            this.coordinates[index+1]['y'],
                                            this.coordinates[index+1]['z'],
                                            pointsPerStep);

    // extend to include these interpolated points, do not include the last
    // element of the array to avoid repeating the number per interval
    interpolatedCoordinatesBuffer = interpolatedCoordinatesBuffer.concat(
        currInterpolation.slice(0, -1));

    // extend the interval buffer
    // credit goes to http://stackoverflow.com/a/13735425/379593
    intervalBuffer = intervalBuffer.concat(
            Array.apply(null, new Array(pointsPerStep)).map(
                Number.prototype.valueOf, index));

  }

  // add the last point to make sure the trajectory is closed
  this.interpolatedCoordinates = interpolatedCoordinatesBuffer.concat(
          currInterpolation.slice(-1));
  this._intervalValues = intervalBuffer;

  return;
}

/**
 *
 * Helper method to calculate the number of points that there should be for a
 * differential.
 *
 * @param {delta} float value for which to determine the required number of
 * points.
 * @param {suppliedN} int value, usually a constant defined by the
 * TrajectoryOfSamples class.
 * @param {minimumDelta} float value to represent the minimum differential
 * found in a set of trajectories.
 *
 * @return an integer representing the number of suggested frames for the
 * differential
 *
 */
TrajectoryOfSamples.prototype.calculateNumberOfPointsForDelta = function(delta){
  return Math.floor((delta*this.suppliedN)/this.minimumDelta);
}

/**
 *
 * Retrieve the representative coordinates needed for a trajectory to be drawn.
 *
 * @param {idx} int value for which to determine the required number of
 * points.
 *
 * @return an Array containing the representative coordinates needed to draw
 * a trajectory at the given index.
 *
 * Note that this implementation is naive and will return points that lay on a
 * rect line if these were part of the original set of coordinates.
 *
 */
TrajectoryOfSamples.prototype.representativeCoordinatesAtIndex = function(idx){

  if (idx === 0){
    return [this.coordinates[0]];
  }

  // we only need to show the edges and none of the interpolated points
  if (this.interpolatedCoordinates.length-1 <= idx){
    return this.coordinates;
  }

  var output = this.coordinates.slice(0, this._intervalValues[idx]+1);
  output = output.concat(this.interpolatedCoordinates[idx]);

  return output;
}

/**
 * 
 * Function to interpolate a certain number of steps between two three
 * dimensional points.
 *
 * This code is based on the function found in:
 *     http://snipplr.com/view.php?codeview&id=47206
 *
 * @param x_1 float initial value of a position in the first dimension
 * @param y_1 float initial value of a position in the second dimension
 * @param z_1 float initial value of a position in the third dimension
 * @param x_2 float final value of a position in the first dimension
 * @param y_2 float final value of a position in the second dimension
 * @param z_2 float final value of a position in the third dimension
 * @param steps integer number of steps that we want the interpolation to run
 *
 * @return Array with a objects that have an x, y and z attributes
 *
 */

function linearInterpolation( x_1, y_1, z_1, x_2, y_2, z_2, steps){
  var xAbs = Math.abs(x_1-x_2);
  var yAbs = Math.abs(y_1-y_2);
  var zAbs = Math.abs(z_1-z_2);
  var xDiff = x_2-x_1;
  var yDiff = y_2-y_1;
  var zDiff = z_2-z_1;

  // and apparetnly this makes takes no effect whatsoever
  var length = Math.sqrt(xAbs*xAbs + yAbs*yAbs + zAbs*zAbs);
  var xStep = xDiff/steps;
  var yStep = yDiff/steps;
  var zStep = zDiff/steps;

  var newx = 0;
  var newy = 0;
  var newz = 0;
  var result = new Array();

  for( var s = 0; s <= steps; s++ ){
    newx = x_1+(xStep*s);
    newy = y_1+(yStep*s);
    newz = z_1+(zStep*s);

    result.push({'x': newx, 'y': newy, 'z': newz});
  }

  return result;
}

/**
 *
 * Function to compute the distance between two three dimensional points.
 *
 * This code is based on the function found in:
 *     http://snipplr.com/view.php?codeview&id=47207
 *
 * @param x_1 float initial value of a position in the first dimension
 * @param y_1 float initial value of a position in the second dimension
 * @param z_1 float initial value of a position in the third dimension
 * @param x_2 float final value of a position in the first dimension
 * @param y_2 float final value of a position in the second dimension
 * @param z_2 float final value of a position in the third dimension
 *
 * @return floating point value of the distance between the two points
 *
 */
function distanceBetweenPoints( x_1, y_1, z_1, x_2, y_2, z_2){
  var xs = 0;
  var ys = 0;
  var zs = 0;

  // Math.pow is faster than simple multiplication
  xs = Math.pow(Math.abs(x_2-x_1), 2);
  ys = Math.pow(Math.abs(y_2-y_1), 2);
  zs = Math.pow(Math.abs(z_2-z_1), 2);

  return Math.sqrt(xs+ys+zs);
}

/**
 *
 * Helper data wrangling function, takes as inputs a mapping file and the
 * coordinates to synthesize the information into an array. Mainly used by the
 * AnimationDirector object.
 * 
 * @param {mappingFileHeaders} an Array of strings containing metadata mapping
 * file headers (required).
 * @param {mappingFileData} an Array where the indices are sample identifiers
 * and each of the contained elements is an Array of strings where the first
 * element corresponds to the first data for the first column in the mapping
 * file (mappingFileHeaders) (required).
 * @param {coordinatesData} an Array of Objects where the indices are the
 * sample identifiers and each of the objects has the following properties: x,
 * y, z, name, color, P1, P2, P3, ... PN where N is the number of dimensions in
 * this dataset (required).
 * @param {trajectoryCategory} a string with the name of the mapping file
 * header where the data that groups the samples is contained, this will
 * usually be BODY_SITE, HOST_SUBJECT_ID, etc. (required).
 * @param {gradientCategory} a string with the name of the mapping file header
 * where the data that spreads the samples over a gradient is contained,
 * usually time or days_since_epoch. Note that this should be an all numeric
 * category (required). 
 *
 * @return an array with the contained data indexed by the sample identifiers.
 *
 */
function getSampleNamesAndDataForSortedTrajectories(mappingFileHeaders,
                                                    mappingFileData,
                                                    coordinatesData,
                                                    trajectoryCategory,
                                                    gradientCategory){
  var gradientIndex = mappingFileHeaders.indexOf(gradientCategory);
  var trajectoryIndex = mappingFileHeaders.indexOf(trajectoryCategory);

  var chewedSampleData = new Object();
  var trajectoryBuffer = null, gradientBuffer = null;

  // this is the most utterly annoying thing ever
  if (gradientIndex === -1) {
    throw new Error("Gradient category not found in mapping file header");
  }
  if (trajectoryIndex === -1) {
    throw new Error("Trajectory category not found in mapping file header");
  }

  for (var sampleId in mappingFileData){

    trajectoryBuffer = mappingFileData[sampleId][trajectoryIndex];
    gradientBuffer = mappingFileData[sampleId][gradientIndex];

    // check if there's already an element for this trajectory, if not
    // initialize a new array for this element of the processed data
    if (chewedSampleData[trajectoryBuffer] === undefined){
      chewedSampleData[trajectoryBuffer] = new Array();
    }
    chewedSampleData[trajectoryBuffer].push({'name': sampleId,
        'value': gradientBuffer, 'x': coordinatesData[sampleId]['x'],
        'y': coordinatesData[sampleId]['y'],
        'z': coordinatesData[sampleId]['z']});
  }

  // we need this custom sorting function to make the values be sorted in
  // ascending order but accounting for the data structure that we just built
  var sortingFunction = function (a, b){
    return parseFloat(a["value"]) - parseFloat(b["value"]);
  }

  // sort all the values using the custom anonymous function
  for (var key in chewedSampleData){
    chewedSampleData[key].sort(sortingFunction);
  }

  return chewedSampleData;
}

/**
 *
 * Function to calculate the minimum delta from an array of wrangled data by
 * getSampleNamesAndDataForSortedTrajectories.
 *
 * @param sampleData is an Array as computed from mapping file data and
 * coordinates by getSampleNamesAndDataForSortedTrajectories.
 *
 * @return float value with the minimum difference between two samples across
 * the defined gradient.
 * 
 * This function will raise an Error in case the input data is undefined.
 * This function will not take into account as a minimum delta zero values i.
 * e.  the differential between two samples that lie at the same position in
 * the gradient.
 *
 */
function getMinimumDelta(sampleData){
  if (sampleData === undefined){
    throw new Error("The sample data cannot be undefined");
  }

  var bufferArray = new Array(), deltasArray = new Array();

  // go over all the data and compute the deltas for all trajectories
  for (var key in sampleData){
    for (var index = 0; index < sampleData[key].length; index++){
      bufferArray.push(sampleData[key][index]['value']);
    }
    for (var index = 0; index < bufferArray.length-1; index++){
      deltasArray.push(Math.abs(bufferArray[index+1]-bufferArray[index]));
    }

    // clean buffer array
    bufferArray.length = 0;
  }

  // remove all the deltas of zero so we don't skew our results
  deltasArray = _.filter(deltasArray, function(num){ return num !== 0; });

  // return the minimum of these values
  return _.min(deltasArray);
}
