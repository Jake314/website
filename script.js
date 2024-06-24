window.addEventListener("load", init)

function init() {
    /* const slide_view_scale = 0.8; */
    var s = 0;
    function scrollFunction(s) {
        $("#scrollingBox").animate({
            /* scrollLeft: $(window).width() * slide_view_scale * s */
            scrollLeft: document.getElementById("scrollingBox").getBoundingClientRect().width * s
        }, 2000);
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