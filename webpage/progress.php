<?php

$cwd[__FILE__] = __FILE__;
if (is_link($cwd[__FILE__])) $cwd[__FILE__] = readlink($cwd[__FILE__]);
$cwd[__FILE__] = dirname($cwd[__FILE__]);

require_once($cwd[__FILE__] . "/../../citizen_science_grid/header.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/navbar.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/footer.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/my_query.php");

$sampler_id = mysql_real_escape_string($_GET['id']);

print_header("DNA@Home: Search $sampler_id Progress", "", "dna");
print_navbar("Projects: DNA@Home", "DNA@Home", "..");

echo "
    <div class='container'>
        <div class='row'>
            <div class='col-sm-12'>";

$sampler_result = query_boinc_db("SELECT samples FROM gibbs_sampler WHERE id = $sampler_id");
$sampler_row = $sampler_result->fetch_assoc();

if (!$sampler_row) {
    echo "Gibbs Sampler run with id $sampler_id not in database.\n";

} else {
    $max_samples = $sampler_row['samples'];

    $walk_result = query_boinc_db("SELECT id, current_steps, had_error FROM gibbs_walk WHERE sampler_id = $sampler_id ORDER BY current_steps DESC");

    while ($walk_row = $walk_result->fetch_assoc()) {
        if ($max_samples == 0) $max_samples = $walk_row['current_steps'];

        $walk_row['progress_percentage'] = ($walk_row['current_steps'] / $max_samples) * 100.0;
//        echo "progress_percentage: " . $walk_row['progress_percentage'] . "<br>";

        if (!$walk_row['had_error']) unset($walk_row['had_error']);

        $walk_rows['row'][] = $walk_row;
    }

    $projects_template = file_get_contents($cwd[__FILE__] . "/templates/walk_overview.html");

    $m = new Mustache_Engine;
    echo $m->render($projects_template, $walk_rows);
}


echo "
            </div> <!-- col-sm-12 -->
        </div> <!-- row -->
    </div> <!-- /container -->";


print_footer('Travis Desell and the DNA@Home Team', 'Travis Desell, Archana Dhasarathy, Sergei Nechaev');

echo "</body></html>";

?>
