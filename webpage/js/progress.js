$(document).ready(function () {

        $('#toggle-active-button:not(.click-bound)').addClass('click-bound').click(function(ev) {
            if ($(this).hasClass('active')) {
                $(this).removeClass('active');
                $(this).text('Hide Active');
                $(".active-progress").removeClass('hidden');
            } else {
                $(this).addClass('active');
                $(this).text('Show Active');
                $(".active-progress").addClass('hidden');
            }
        }); 

        $('#toggle-finished-button:not(.click-bound)').addClass('click-bound').click(function(ev) {
            if ($(this).hasClass('active')) {
                $(this).removeClass('active');
                $(this).text('Hide Finished');
                $(".finished-progress").removeClass('hidden');
            } else {
                $(this).addClass('active');
                $(this).text('Show Finished');
                $(".finished-progress").addClass('hidden');
            }
        }); 

        $('#toggle-error-button:not(.click-bound)').addClass('click-bound').click(function(ev) {
            if ($(this).hasClass('active')) {
                $(this).removeClass('active');
                $(this).text('Hide Errors');
                $(".error-progress").removeClass('hidden');
            } else {
                $(this).addClass('active');
                $(this).text('Show Errors');
                $(".error-progress").addClass('hidden');
            }
        }); 


});

