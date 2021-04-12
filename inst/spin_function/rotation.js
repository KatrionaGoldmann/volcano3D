/* Toggles spinning based on action button event
*/
function plot_rotate(el) {
  // Update the colour of the rotate button
  $("[data-title='Rotate'] svg path").css("fill", "#7ac143");
  
   // get the last button added (-2 because of plotly button and 0 index)
  let buttons = document.getElementsByClassName("modebar-btn");
  let titles = Array();
  for (let i = 0; i < buttons.length; i++) {
     titles.push(buttons[i].getAttribute('data-title'));
  }
  var useIndex = titles.indexOf('Rotate');
  
  //Rotate the plot  
  function run() {
    var modeBarButtons = document.getElementsByClassName("modebar-btn")[useIndex];
    var shouldRotate = (modeBarButtons.getAttribute("data-val") == 'true');

    if (shouldRotate) {
      rotate("scene", Math.PI / 180);
    }
    requestAnimationFrame(run);  // Recursive call
  }

  run(); // Start the animation

  function rotate(id, angle) {
    var eye0 = el.layout[id].camera.eye;
    console.log(eye0);
    if(! eye0) {
      eye0 = {
        x: 1,
        y: 1,
        z: 1
      };
    }
    console.log(eye0);
    
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
