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
 * @name THREE.EmperorTrajectory
 *
 * @class This class represents the internal logic for a linearly interpolated
 * tube/trajectory in THREE.js the object itself is a subclass of the
 * THREE.Curve.
 *
 * @credits: This answer in StackOverflow helped a lot:
 * http://stackoverflow.com/a/18580832/379593
 *
 */
THREE.EmperorTrajectory = THREE.Curve.create(
  function ( points) {
    this.points = (points == undefined) ? [] : points;
  },

  function ( t ) {    
    var points = this.points;
    var index = ( points.length - 1 ) * t;
    var floorIndex = Math.floor(index);

    if(floorIndex == points.length-1){
      return points[floorIndex];
    }

    var floorPoint = points[floorIndex];
    var ceilPoint = points[floorIndex+1];

    return floorPoint.clone().lerp(ceilPoint, index - floorIndex);
  }
);

THREE.EmperorTrajectory.prototype.getUtoTmapping = function(u) {
    return u;
};

/**
 *
 * Format an SVG string with labels and colors.
 *
 * @param {labels} Array object with the name of the labels.
 * @param {colors} Array object with the colors for each label.
 *
 * @return Returns an SVG string with the labels and colors values formated as
 * a legend.
 *
 */
function formatSVGLegend(labels, colors){
  var labels_svg='', pos_y=1, increment=40, max_len=0, rect_width,
      font_size=12;

  for (var i=0; i<labels.length; i++){
    // add the rectangle with the corresponding color
    labels_svg += '<rect height="27" width="27" y="'+pos_y+
        '" x="5" style="stroke-width:1;stroke:rgb(0,0,0)" fill="'+
        colors[i]+'"/>';

    // add the name of the category
    labels_svg += '<text xml:space="preserve" y="'+(pos_y+20)+'" x="40" '+
        'font-size="'+font_size+'" stroke-width="0" stroke="#000000" '+
        'fill="#000000">'+labels[i]+'</text>';

    pos_y += increment;
  }

  // get the name with the maximum number of characters and get the length
  max_len = _.max(labels, function(a){return a.length}).length

  // duplicate the size of the rectangle to make sure it fits the labels
  rect_width = font_size*max_len*2;

  labels_svg = '<svg xmlns="http://www.w3.org/2000/svg" width="'+
      rect_width+'" height="'+(pos_y-10)+'"><g>'+labels_svg+'</g></svg>';

  return labels_svg;
}
