<?php

$cwd[__FILE__] = __FILE__;
if (is_link($cwd[__FILE__])) $cwd[__FILE__] = readlink($cwd[__FILE__]);
$cwd[__FILE__] = dirname($cwd[__FILE__]);

require_once($cwd[__FILE__] . "/../../citizen_science_grid/header.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/navbar.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/footer.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/my_query.php");

print_header("DNA@Home: Current Search Progress", "<script type='text/javascript' src='./js/overview.js'></script>", "dna");
print_navbar("Projects: DNA@Home", "DNA@Home", "..");

echo "
    <div class='container'>
        <div class='row' style='margin-bottom:10px;'>
            <div class='col-sm-12'>
                <button type='button' id='display-inactive-runs-button' class='btn btn-primary pull-right'>
                    Display Inactive Runs
                </button>
            </div> <!--col-sm-12 -->
        </div> <!-- row -->

        <div class='row'>
            <div class='col-sm-12'>";

$sampler_result = query_boinc_db("SELECT * FROM gibbs_sampler ORDER BY name");

while ($sampler_row = $sampler_result->fetch_assoc()) {
    if ($sampler_row['samples'] == 0) $sampler_row['is_hidden'] = true;

    $walks_result = query_boinc_db("SELECT AVG(current_steps) FROM gibbs_walk WHERE sampler_id = " . $sampler_row['id']);
    $walks_row = $walks_result->fetch_assoc();

    $sampler_row['avg_progress'] = $walks_row['AVG(current_steps)'];

    $sampler_rows['row'][] = $sampler_row;
}

$projects_template = file_get_contents($cwd[__FILE__] . "/templates/sampler_overview.html");

$m = new Mustache_Engine;
echo $m->render($projects_template, $sampler_rows);


echo "
            </div> <!-- col-sm-12 -->
        </div> <!-- row -->

    </div> <!-- /container -->";


print_footer('Travis Desell and the DNA@Home Team', 'Travis Desell, Archana Dhasarathy, Sergei Nechaev');

echo "</body></html>";

?>
