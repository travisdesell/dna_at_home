$(document).ready(function () {

        $('#display-inactive-runs-button:not(.click-bound)').addClass('click-bound').click(function(ev) {
            if ($(this).hasClass('active')) {
                $(this).removeClass('active');
                $(this).text('Display Inactive Runs');
            } else {
                $(this).addClass('active');
                $(this).text('Hide Inactive Runs');
            }

            $(".inactive-search").toggle();
        }); 

});

