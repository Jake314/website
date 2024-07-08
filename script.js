window.addEventListener("load", init)

function init() {
    var s = 0;
    function scrollFunction(s) {
        $("#scrollingBox").animate({
            scrollLeft: document.getElementById("scrollingBox").getBoundingClientRect().width * s
        }, 1000);
    };

    $("#rightButton").click(function() {
        s = (s + 1) % 6;
        scrollFunction(s);
      });
    $("#leftButton").click(function() {
        s -= 1;
        if (s < 0) {
            s += 6;
        }
        scrollFunction(s);
    });
}