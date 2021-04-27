/* Toggles spinning based on action button event
*/
function plot_rotate(el) {
  // Update the colour of the rotate button
  $("[data-title='Rotate'] svg path").css("fill", "rotate_colour");
  
   // get the index for the rotate button
  let buttons = document.getElementsByClassName("modebar-btn");
  let titles = Array();
  for (let i = 0; i < buttons.length; i++) {
     titles.push(buttons[i].getAttribute('data-title'));
  }
  var useIndex = titles.indexOf('Rotate');
  
  function run() {
    rotate("scene", Math.PI / rotation_speed);
    spinID = requestAnimationFrame(run);  // Recursive call
  }

  // Start/stop the animation
  $("[data-title='Rotate']").on("click", function() {
  	var modeBarButtons = document.getElementsByClassName("modebar-btn")[useIndex];
    var shouldRotate = (modeBarButtons.getAttribute("data-val") == 'true');
    if (shouldRotate) {
  	  spinID = requestAnimationFrame(run);
  	} else {
      cancelAnimationFrame(spinID);
    }
  });
  
  // stop animation if change in plot type
  $(document).on('shiny:inputchanged', function(event) {
    if([shiny_event_names].includes(event.name)) {
      cancelAnimationFrame(spinID);
    }
  });

  function rotate(id, angle) {
    var eye0 = el.layout[id].camera.eye;
    if(! eye0) {
      eye0 = {
        x: 1,
        y: 1,
        z: 1
      };
    }
    
    var rtz = xyz2rtz(eye0);
    rtz.t += angle;
    var eye1 = rtz2xyz(rtz);
    Plotly.relayout(el, id + ".camera.eye", eye1);
  }

  function xyz2rtz(xyz) {
    return {
      r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
      t: Math.atan2(xyz.y, xyz.x),
      z: xyz.z,
    };
  }

  function rtz2xyz(rtz) {
    return {
      x: rtz.r * Math.cos(rtz.t),
      y: rtz.r * Math.sin(rtz.t),
      z: rtz.z,
    };
  }
}
