var vOffset_init = 65;
var vOffset = vOffset_init;
var c = 'collapsed';

function toggleList(toggle, content, maxItems) {
    if (toggle.css('display') == 'none') {
        vOffset = vOffset_init;
        toggle.removeClass(c);
        content.show();
        return;
    } else
        vOffset = 8;

    if (maxItems > content.children().length)
        return;
    content.hide();
    toggle.addClass(c);
}

$(function () {
    $('a[href*=#]:not([href=#])').on('click', function (e) {
        if (e.which == 2 || e.metaKey || e.ctrlKey || e.shiftKey)
            return true;
        var target = $(this.hash.replace(/(\.)/g, "\\$1"));
        target = target.length ? target : $('[name=' + this.hash.slice(1) + ']');
        if (target.length) {
            setTimeout(function () {
                $('html, body').animate({scrollTop: target.offset().top - vOffset}, 50);}, 50);
        }
    });
});

$(window).load(function () {
    var hashChanged = function() {
        var h = window.location.hash;
        var re = /[^a-z0-9_\.\#\-]/i
        if (h.length > 1 && !re.test(h)) {
            setTimeout(function () {
                var tgt = $(h.replace(/(\.)/g, "\\$1"));
                tgt = tgt.length ? tgt : $('[name=' + h.slice(1) + ']');
                $(window).scrollTop(tgt.offset().top - vOffset);
            }, 0);
        }
    }
    $(window).bind('hashchange', hashChanged);
    hashChanged.call();

    if (!$('.sidebar toc').is(':empty')) {
        $('<div id="toc-toggle"></div>').prependTo('.sidebar .toc');
        var toc = $('.sidebar .toc ul');
        var tocToggle = $('#toc-toggle');
        var tocCallback = function() { toggleList(tocToggle, toc, 4); };

        $('#toc-toggle').on('click', function(e) {
            e.stopPropagation();
            toc.toggle();
            tocToggle.toggleClass(c);
        });

        tocCallback.call();
        $(window).resize(tocCallback);
    }

    if (!$('#sidebar-content').is(':empty')) {
        $('#sidebar-content h2').first().clone().prependTo('#sidebar-content');
        $('<div id="sidebar-toggle"></div>').prependTo('#sidebar-content');
        var sb = $('#sidebar-content .sectionlist');
        var sbToggle = $('#sidebar-toggle');
        var sbCallback = function() { toggleList(sbToggle, sb, 0); };

        $('#sidebar-toggle').on('click', function(e) {
            e.stopPropagation();
            sb.toggle();
            sbToggle.toggleClass(c);
        });

        sbCallback.call();
        $(window).resize(sbCallback);
    }
});
