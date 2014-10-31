<?php

$cwd[__FILE__] = __FILE__;
if (is_link($cwd[__FILE__])) $cwd[__FILE__] = readlink($cwd[__FILE__]);
$cwd[__FILE__] = dirname($cwd[__FILE__]);

require_once($cwd[__FILE__] . "/../../citizen_science_grid/header.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/navbar.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/footer.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/my_query.php");

$sampler_id = mysql_real_escape_string($_GET['id']);

$sampler_result = query_boinc_db("SELECT samples, name FROM gibbs_sampler WHERE id = $sampler_id");
$sampler_row = $sampler_result->fetch_assoc();

if (!$sampler_row) {
    echo "Gibbs Sampler run with id $sampler_id not in database.\n";

    print_header("DNA@Home: Progress for Unknown Search", "", "dna");
    print_navbar("Projects: DNA@Home", "DNA@Home", "..");
    print_footer('Travis Desell and the DNA@Home Team', 'Travis Desell, Archana Dhasarathy, Sergei Nechaev');


} else {
    $sampler_name = $sampler_row['name'];

    print_header("DNA@Home: Progress for $sampler_name", "<script type='text/javascript' src='./js/progress.js'></script>", "dna");
    print_navbar("Projects: DNA@Home", "DNA@Home", "..");

    echo "
        <div class='container'>
        <div class='row' style='margin-bottom:10px;'>
            <div class='col-sm-12'>
                <a type='button' class='btn btn-default pull-left' href='./overview.php'>
                    Return to Overview
                </a>

                <button id='toggle-finished-button' type='button' class='btn btn-success pull-right' href='./overview.php'>
                    Hide Finished
                </a>
                <button id='toggle-active-button' type='button' class='btn btn-primary pull-right' href='./overview.php'>
                    Hide Active
                </a>
                <button id='toggle-error-button' type='button' class='btn btn-danger pull-right' href='./overview.php'>
                    Hide Errors
                </button>
            </div> <!-- col-sm-12 -->
        </div> <!-- row -->

        <div class='row'>
        <div class='col-sm-12'>";

    $max_samples = $sampler_row['samples'];

    $walk_result = query_boinc_db("SELECT id, current_steps, had_error FROM gibbs_walk WHERE sampler_id = $sampler_id ORDER BY current_steps DESC");

    while ($walk_row = $walk_result->fetch_assoc()) {
        if ($max_samples == 0) $max_samples = $walk_row['current_steps'];

        $walk_row['progress_percentage'] = ($walk_row['current_steps'] / $max_samples) * 100.0;
//        echo "progress_percentage: " . $walk_row['progress_percentage'] . "<br>";

        if ($walk_row['progress_percentage'] == 100) $walk_row['finished'] = true;
        else $walk_row['active'] = true;

        $walk_row['max_samples'] = $max_samples;
        if (!$walk_row['had_error']) unset($walk_row['had_error']);
        else unset($walk_row['active']);

        $walk_rows['row'][] = $walk_row;
    }

    $projects_template = file_get_contents($cwd[__FILE__] . "/templates/walk_overview.html");

    $m = new Mustache_Engine;
    echo $m->render($projects_template, $walk_rows);


    echo "
        </div> <!-- col-sm-12 -->
        </div> <!-- row -->
        </div> <!-- /container -->";


    print_footer('Travis Desell and the DNA@Home Team', 'Travis Desell, Archana Dhasarathy, Sergei Nechaev');

}

echo "</body></html>";

?>
