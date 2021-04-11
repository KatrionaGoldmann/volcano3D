/* Toggles spinning based on action button event
*/
function(el) {
  $("[data-title='Rotate'] svg path").css("fill", "#7ac143");
  function run() {
    var modeBarButtons = document.getElementsByClassName("modebar-btn")[7];
    var shouldRotate = eval(modeBarButtons.getAttribute("data-val"));

    if (shouldRotate) {
      rotate("scene", Math.PI / 180);
    }
    requestAnimationFrame(run);  // Recursive call
  }

  run(); // Start the animation

  function rotate(id, angle) {
    var eye0 = el.layout[id].camera.eye;
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
