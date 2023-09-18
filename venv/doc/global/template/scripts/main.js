"use strict";

function createCookie(name, value, days) {
    var expires;
    if (days) {
        var date = new Date();
        date.setTime(date.getTime() + (days * 24 * 60 * 60 * 1000));
        expires = "; expires=" + date.toGMTString();
    } else {
        expires = "";
    }
    document.cookie = escape(name) + "=" + escape(value) + expires + "; path=/";
    $('.cookies_yum').click(function() {
        $(this).fadeOut()
    });
}
function readCookie(name) {
    var nameEQ = escape(name) + "=";
    var ca = document.cookie.split(';');
    for (var i = 0; i < ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0) === ' ') c = c.substring(1, c.length);
        if (c.indexOf(nameEQ) === 0) return unescape(c.substring(nameEQ.length, c.length));
    }
    return null;
}
function eraseCookie(name) {
    createCookie(name, "", -1);
}
function load_sdk(s, id, src) {
    var js, fjs = document.getElementsByTagName(s)[0];
    if (document.getElementById(id)) return;
    js = document.createElement(s);
    js.id = id;
    js.src = src;
    fjs.parentNode.insertBefore(js, fjs);
}
$(document).ready(function($) {
    if (document.documentElement.clientWidth < 1280) {
        oneQt.extraLinksToMain();
    }

    $('#menuextras .search').click(function(e){
        e.preventDefault();
        $('.big_bar.account').slideUp();
        $('.big_bar.search').slideToggle();
        $('.big_bar_search').focus();
        $(this).toggleClass('open');
    });
    $('.cookies_yum').click(function() {
        $('.cookies_yum').fadeOut();
        createCookie("cookies_nom", "yum", 180);
        var cookie_added = 1;
    });
    if (!(readCookie('cookies_nom') == 'yum')) {
        $('.cookies_yum').fadeIn();
    } else {
        var cookie_added = 1;
    }

    Modernizr.load({test: Modernizr.input.placeholder,
                    nope: wpThemeFolder + '/js/placeholders.min.js'});

    $('#navbar .navbar-toggle').click(function(e) {
        e.preventDefault();
        if ($(this).hasClass('opened')) {
            $(this).removeClass('opened');
            $('#navbar .navbar-menu').css('max-height', '0px');
        }
        else {
            $(this).addClass('opened');
            $('#navbar .navbar-menu').css('max-height', $('#navbar .navbar-menu ul').outerHeight() + 'px');
        }
    });

    $(window).resize(function() {
        oneQt.stickySidebar();
        oneQt.footerPosition();
        if (document.documentElement.clientWidth < 1280) {
            oneQt.extraLinksToMain();
        } else {
            oneQt.mainLinkstoExtra();
        }
    });

    $(window).scroll(function() {
        oneQt.stickySidebar();
        oneQt.stickyHeader();
    });

    oneQt.stickySidebar();
    oneQt.footerPosition();
    oneQt.tabContents();
});

$( window ).load(function() {
    load_sdk('script', 'facebook-jssdk','//connect.facebook.net/en_US/sdk.js#xfbml=1&appId=207346529386114&version=v2.0');
    load_sdk('script', 'twitter-wjs', '//platform.twitter.com/widgets.js');
    $.getScript("//www.google.com/jsapi", function(){
        google.load("feeds", "1", {"callback": oneQt.liveFeeds});
    });
});

var oneQt = {
    stickySidebar: function() {
        if ($('#sidebar').length && $('#sidebar').outerHeight() > 20) {
            var $sidebar = $('#sidebar');
            var $win = $(window);
            var $sidebarContainer = $sidebar.parent();
            var headerHeight = $('#navbar').outerHeight();
            if ($win.outerHeight() - headerHeight > $sidebar.innerHeight() &&
                    $win.scrollTop() > $sidebarContainer.offset().top) {
                var newTop = headerHeight + $win.scrollTop() - $sidebarContainer.offset().top;
                if (newTop + $sidebar.innerHeight() > $sidebarContainer.innerHeight())
                    newTop = $sidebarContainer.innerHeight() - $sidebar.innerHeight();

                $sidebar.css({top: newTop +'px'})
            }
            else {
                $sidebar.css({top: '0'})
            }
        }
    },

    footerPosition: function () {
        $('#footerbar').removeClass('fixed');
        if (($('.hbspt-form').length > 0) || ($('#customerInfo').length > 0) || ($('.purchase_bar').length > 0)) {
            var footerBottomPos = $('#footerbar').offset().top + $('#footerbar').outerHeight();
            if (footerBottomPos < $(window).height())
                $('#footerbar').addClass('fixed');
        }
    },

    stickyHeader: function () {
        var originalHeaderHeight = 79;
        if ($(window).scrollTop() > originalHeaderHeight) {
            $('#navbar').addClass('fixed');
            $('#bottom_header').fadeOut();

            if (!(cookie_added == 1)) {
                $('.cookies_yum').fadeOut();
                createCookie("cookies_nom", "yum", 180);
                var cookie_added = 1;
            }
        }
        else {
            $('#navbar').removeClass('fixed');
            $('#bottom_header').fadeIn();
        }
    },

    tabContents: function () {
        $('.tab-container').each(function(i) {
            var $el = $(this);
            $el.find('.tab-titles li:eq(0)').addClass('active');
            $el.find('.tab-contents .tab:eq(0)').addClass('active');
            $el.find('.tab-titles a').click(function(e) {
                e.preventDefault();
                var index = $(this).parent().index();
                $el.find('.tab-titles li').removeClass('active');
                $el.find('.tab-contents .tab').removeClass('active');
                $(this).parent().addClass('active');
                $el.find('.tab-contents .tab').eq(index).addClass('active');
            })
        });
    },

    liveFeeds: function () {
        $('.feed-container').each(function(i) {
            var feedUrl = $(this).data('url');
            if (feedUrl != "") oneQt.blogFeed($(this), feedUrl);
        });
    },

    blogFeed: function ($container, feedUrl) {
        var feed = new google.feeds.Feed(feedUrl);
        feed.setNumEntries(3);
        feed.load(function(result) {
            $container.html('');
            if (!result.error) {
                for (var i = 0; i < result.feed.entries.length; i++) {
                    var entry = result.feed.entries[i];
                    var $article = $('<article class="discussion-tile cf"></article>');
                    $container.append($article);
                    var html = '    <div class="author retina">';
                    html += '        <img src="'+wpThemeFolder+'/assets/images/author_placeholder.png" alt="">';
                    html += '    </div>';
                    html += '    <div class="discussion-item">';
                    html += '        <h4><a href="'+encodeURI(entry.link)+'"></a></h4>'
                    html += '        <h3><a href="'+encodeURI(entry.link)+'" target="_blank"></a></h3>'
                    html += '        <p><a href="'+encodeURI(entry.link)+'" target="_blank"></a></p>';
                    html += '        <ul class="taglist cf">';
                    html += '        </ul>';
                    html += '    </div>';
                    $article.append(html);
                    $article.find('h4 a').text(result.feed.title);
                    $article.find('h3 a').text(entry.title);
                    $article.find('p a').text(entry.author);
                    try {
                        for (var j=0; j<entry.categories.length; j++) {
                            var $li = $('<li><a href="'+encodeURI(entry.link)+'" target="_blank" class="btn btn-tag"></a></li>');
                            $li.find('a').text(entry.categories[j]);
                            $article.find('.taglist').append($li);
                        }
                    } catch(e) {}
                }
                if (result.feed.link && result.feed.link != "") {
                    var linkHtml = '<a href="'+encodeURI(result.feed.link)+'" class="text-lightgrey" target="_blank">Show all</a>';
                    $container.append(linkHtml);
                }
            }
        });
    },

    extraLinksToMain: function() {
        var extramenuLinks = $('#menuextras').find('li');
        var mainmenu = $('#mainmenu');
        var count = 0;
        if ($(extramenuLinks).length > 2) {
            $(extramenuLinks).each(function() {
                if (count < 3) {
                    var newLink = $(this);
                    $(newLink).addClass('dynamic-add');
                    $(mainmenu).append(newLink);
                }
                count++;
            });
        }
    },

    mainLinkstoExtra: function() {
        var mainmenuLinks = $('#mainmenu').find('.dynamic-add');
        var extramenu = $('#menuextras');
        var count = 0;
        $(mainmenuLinks).each(function() {
            var newLink = $(this);
            $(extramenu).prepend(newLink);
            count++;
        });
    }
}
