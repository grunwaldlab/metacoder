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


// http://colorbrewer2.org > qualitative > Number of Data Classes = 12
// colorbrewer will provide you with two lists of colors, those have been
// added here
var colorbrewerDiscrete = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
    "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
    "#ccebc5", "#ffed6f", /*first list ends here*/
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"];
// taken from the qiime/colors.py module; a total of 24 colors
var qiimeDiscrete = [ "#ff0000", "#0000ff", "#f27304", "#008000", "#91278d",
    "#ffff00", "#7cecf4", "#f49ac2", "#5da09e", "#6b440b", "#808080",
    "#f79679", "#7da9d8", "#fcc688", "#80c99b", "#a287bf", "#fff899",
    "#c49c6b", "#c0c0c0", "#ed008a", "#00b6ff", "#a54700", "#808000",
    "#008080"];
var k_DiscreteColorMaps = {"discrete-coloring":colorbrewerDiscrete,
                           "discrete-coloring-qiime":qiimeDiscrete};

// these colors are included in chroma and are the only ones we should
// use to interpolate through whe coloring in a continuous mode
var k_CHROMABREWER_MAPS = ['discrete-coloring-qiime', 'discrete-coloring',
    'OrRd', 'PuBu', 'BuPu', 'Oranges', 'BuGn', 'YlOrBr', 'YlGn', 'Reds',
    'RdPu', 'Greens', 'YlGnBu', 'Purples', 'GnBu', 'Greys', 'YlOrRd', 'PuRd',
    'Blues', 'PuBuGn', 'Spectral', 'RdYlGn', 'RdBu', 'PiYG', 'PRGn', 'RdYlBu',
    'BrBG', 'RdGy', 'PuOr'];
var k_CHROMABREWER_MAPNAMES = ['Classic QIIME Colors',
    'Discrete Coloring (Colorbrewer)', 'Orange-Red', 'Purple-Blue',
    'Blue-Purple', 'Oranges', 'Blue-Green', 'Yellow-Orange-Brown',
    'Yellow-Green', 'Reds', 'Red-Purple', 'Greens', 'Yellow-Green-Blue',
    'Purples', 'Green-Blue', 'Greys', 'Yellow-Orange-Red', 'Purple-Red',
    'Blues', 'Purple-Blue-Green', 'Spectral', 'Red-Yellow-Green', 'Red-Blue',
    'Pink-Yellow-Green', 'Pink-Red-Green', 'Red-Yellow-Blue',
    'Brown-Blue-Green', 'Red-Grey', 'Purple-Orange'];

/**
 *
 * Sorting function that deals with alpha and numeric elements
 * 
 * This function takes a list of strings, divides it into two new lists, one
 * that's alpha-only and one that's numeric only.
 *
 * @param {list} an Array of strings.
 *
 * @return a sorted Array where all alpha elements at the beginning & all
 * numeric elements at the end.
 *
 */
function naturalSort(list){
  var numericPart = [], alphaPart = [], result = [];

  // separate the numeric and the alpha elements of the array
  for(var index = 0; index < list.length; index++){
    if(isNaN(parseFloat(list[index]))){
      alphaPart.push(list[index])
    }
    else{
      numericPart.push(list[index])
    }
  }

  // ignore casing of the strings, taken from:
  // http://stackoverflow.com/a/9645447/379593
  alphaPart.sort(function (a, b) {
    return a.toLowerCase().localeCompare(b.toLowerCase());
  });

  // sort in ascending order
  numericPart.sort(function(a,b){return parseFloat(a)-parseFloat(b)})

  return result.concat(alphaPart, numericPart);
}

/**
 *
 * Utility function to convert an XML DOM documents to a string useful for unit
 * testing
 *
 * @param {node} XML DOM object, usually as created by the document object.
 *
 * @return string representation of the node object.
 *
 * This code is based on this answer http://stackoverflow.com/a/1750890
 *
 */
function convertXMLToString(node) {
  if (typeof(XMLSerializer) !== 'undefined') {
    var serializer = new XMLSerializer();
    return serializer.serializeToString(node);
  }
  else if (node.xml) {
    return node.xml;
  }
}

/**
 *
 * Retrieve a discrete color.
 *
 * @param {index} int, the index of the color to retrieve.
 * @param {map} string, name of the discrete color map to use.
 *
 * @return string representation of the hexadecimal value for a color in the
 * list the QIIME colors or the ColorBrewer discrete colors. If this value
 * value is greater than the number of colors available, the function will just
 * rollover and retrieve the next available color.
 *
 * Defaults to use ColorBrewer colors if there's no map passed in.
 *
 */
function getDiscreteColor(index, map){
  if (map === undefined){
    map = 'discrete-coloring';
  }
  if (_.has(k_DiscreteColorMaps, map) === false){
    throw new Error("Could not find "+map+" as a discrete colormap.")
  }

  var size = k_DiscreteColorMaps[map].length;
  if(index >= size){
    index = index - (Math.floor(index/size)*size)
  }

  return k_DiscreteColorMaps[map][index]
}


/**
 *
 * Generate a list of colors that corresponds to all the samples in the plot
 *
 * @param {values} list of objects to generate a color for, usually a category
 * in a given metadata column.
 * @param {map} name of the color map to use, see k_CHROMABREWER_MAPS.
 *
 *
 * This function will generate a list of coloring values depending on the
 * coloring scheme that the system is currently using (discrete or continuous).
*/
function getColorList(values, map) {
  var colors = {}, numColors = values.length-1, counter=0, interpolator,
      discrete = false;

  if (k_CHROMABREWER_MAPS.indexOf(map) === -1) {
    throw new Error("Could not find "+map+" in the available colormaps");
  }

  if (map === 'discrete-coloring' || map === 'discrete-coloring-qiime'){
    discrete = true;
  }

  // 1 color and continuous coloring should return the first element of the map
  if (numColors === 0 && discrete === false){
    colors[values[0]] = chroma.brewer[map][0];
    return colors;
  }

  if (discrete === false){
    map = chroma.brewer[map];
    interpolator = chroma.interpolate.bezier([map[0], map[3], map[4], map[5],
                                              map[8]]);
  }

  for(var index in values){
    if(discrete){
      // get the next available color
      colors[values[index]] = getDiscreteColor(index, map);
    }
    else{
      colors[values[index]] =  interpolator(counter/numColors).hex();
      counter = counter + 1;
    }
  }

  return colors;
}


/**
 *
 * Escape special characters in a string for use in a regular expression.
 *
 * @param {regex} string to escape for use in a regular expression.
 *
 * @return string with escaped characters for use in a regular expression.
 *
 * Credits go to this SO answer http://stackoverflow.com/a/5306111
 */
function escapeRegularExpression(regex){
    return regex.replace(/[-[\]{}()*+?.,\\^$|#\s]/g, "\\$&");
}

/**
 *
 * Clean a string in HTML formatted strings that get created with the namespace
 * tag in some browsers and not in others. Intended to facilitate testing.
 *
 * @param {htmlString} string to remove namespace from.
 *
 * @return string without namespace.
 *
 */
function cleanHTML(htmlString){
    return htmlString.replace(' xmlns="http://www.w3.org/1999/xhtml"', '')
}

