function spin_button(gd, ev) {
  var button = ev.currentTarget;
  var val = (button.getAttribute("data-val") == 'true');
  
  // update whether to rotate
  button.setAttribute("data-val", ! val);
  
  // update the title and logo
  var title = button.getAttribute("data-title");
  
  if(title == "Rotate"){
    $("[data-title=\'Rotate\'] svg path").css("fill", "stop_colour");
    button.firstElementChild.firstElementChild.setAttribute("d","stop_icon_path");
    button.setAttribute("data-title", "Stop");
    
  } else {
    $("[data-title=\'Stop\'] svg path").css("fill", "rotate_colour");
    button.firstElementChild.firstElementChild.setAttribute("d","rotate_icon_path");
    button.setAttribute("data-title", "Rotate");
  }
}
