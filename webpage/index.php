<?php

$cwd[__FILE__] = __FILE__;
if (is_link($cwd[__FILE__])) $cwd[__FILE__] = readlink($cwd[__FILE__]);
$cwd[__FILE__] = dirname($cwd[__FILE__]);

require_once($cwd[__FILE__] . "/../../citizen_science_grid/header.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/navbar.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/footer.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/my_query.php");

require_once("/projects/csg/html/inc/text_transform.inc");

print_header("DNA@Home");
print_navbar("", "Citizen Science Grid: DNA@Home");

echo "
    <div class='container'>
      <div class='row'>
        <div class='col-sm-12'>";

include $cwd[__FILE__] . "/templates/dna_info.html";

echo "
        </div> <!-- col-sm-12 -->
    </div> <!-- row -->

    <div class='row'>
        <div class='col-sm-6'>
            <div class='well' style='margin-top:20px;'>
                <h3>User of the Day</h3>

            </div>

            <div class='well'>
                <h3><a href='http://volunteer.cs.und.edu/csg/forum_forum.php?id=1'>News</a> <img href='rss_main.php' src='http://volunteer.cs.und.edu/csg/img/rss_icon.gif' alt='RSS'> </h3></p>";


$thread_result = query_boinc_db("SELECT id, title, owner, timestamp FROM thread WHERE forum = 1 AND hidden = 0 ORDER BY id desc LIMIT 5");


while ($thread_row = $thread_result->fetch_assoc()) {
    $post_result = query_boinc_db("SELECT content, timestamp FROM post WHERE thread = " . $thread_row['id'] . " ORDER BY id LIMIT 1");
    $post_row = $post_result->fetch_assoc();

    $owner = csg_get_user_from_id($thread_row['owner']);

    echo "
                <hr class='news-hr'>
                <p><b>" . $thread_row['title'] . "</b></p>
                <p>" . output_transform($post_row['content']) . "</p>
                <p style='text-align:right; margin-bottom:0px'><i>" . $owner['name'] . " on " . date("l, F jS", $thread_row['timestamp']) . "</i><br>
                <a href='http://volunteer.cs.und.edu/csg/forum_thread.php?id=" . $thread_row['id'] . "'>leave a comment</a></p>";

}

echo "
            </div> <!-- well -->
        </div> <!-- col-sm-6 -->

        <div class='col-sm-6'>
            <div class='well' style='margin-top:20px;'>
                <h3>Notice for DNA@Home Users</h3>
                <p>DNA@Home is transitioning to a sub-project of Citizen Science Grid. DNA@Home workunits will be sent out from the Citizen Science Grid project.  You can link your DNA@Home account to your account on Citizen Science Grid so it will continue to gain credit for completing DNA@Home workunits on Citizen Science Grid by visiting the <a href='../csg/link_accounts.php'>link accounts</a> webpage.</p>
            </div>

            <div class='well' style='margin-top:20px;'>
                <h3>Contact Information</h3>
                DNA@Home is currently run by:
                <ul>
                    <li>
                    <a href='http://people.cs.und.edu/~tdesell/'>Travis Desell</a>, Assistant Professor of Computer Science, University of North Dakota
                    </li>
                    <li>
                    <a href='http://www.med.und.edu/biochemistry/faculty-staff/archana-dhasarathy.cfm'>Archana Dhasarathy</a>, Assistant Professor of Anatomy and Cell Biology, University of North Dakota
                    </li>
                    <li>
                    <a href='http://www.med.und.edu/anatomy/contact/sergei-nechaev.cfm'>Sergei Nechaev</a>, Assistant Professor of Anatomy and Cell Biology, University of North Dakota
                    </li>
                    <li>
                    Kris Zarns, Graduate Research Assistant in Computer Science, University of North Dakota
                    </li>
                </ul>

                With previous support from:
                <ul>
                    <li>
                    <a href='http://www.cs.rpi.edu/~magdon/'>Malik Magdon-Ismail</a>, Associate Professor of Computer Science
                    </li>
                    <li>
                    <a href='http://www.rpi.edu/~newbel/'>Lee Newberg</a> Research Associate Professor of Computer Science, and New York State Wadsworth Center Research Scientist 
                    </li>
                    <li>
                    <a href='http://www.cs.rpi.edu/~szymansk/index.php'>Boleslaw Szymanski</a>, Claire and Roland Schmitt Distinguished Professor of Computer Science 
                    </li>
                    <li>
                    Lei Chen, Graduate Research Assistant in Computer Science, Rensselaer Polytechnic Institute
                    </li>
                </ul>
            </div>

            <div class='well' style='margin-top:20px;'>
                <h3>Publications</h3>
                <p>
                    All publications, talks and posters about DNA@Home can be viewed on the Citizen Science Grid's <a href='../csg/publications.php'>publications page</a>.
                </p>
            </div>

            <div class='well' style='margin-top:20px;'>
                <h3>Support</h3>
                    <div class='row'>
                        <div class='col-sm-3'>
                            <p align=center>
                            <a href='http://und.edu'><img class='img-responsive' src='http://volunteer.cs.und.edu/wildlife/images/und_logo.png'></a>
                            </p>
                        </div>
                        <div class='col-sm-9'>
                        <p>
                            DNA@Home has been generously supported by a new faculty SEED grant from UND's Office of Research Development and Compliance, and a SEED Grant from UND's Medical School. The Citizen Science Grid volunteer computing server is hosted by UND's <a href='http://www.aero.und.edu/about/SCC.aspx'>Scientific Computing Center</a>.
                        </p>
                        </div>
                    </div>

                    <div class='row'>
                        <div class='col-sm-3'>
                            <a href=\"http://boinc.berkeley.edu/\" style='display:block; margin-left:auto; margin-right:auto;'><img class='img-responsive' src=\"http://volunteer.cs.und.edu/wildlife/img/pb_boinc.gif\" alt=\"Powered by BOINC\"></a>
                        </div>
                        <div class='col-sm-9'>
                            <p>
                            DNA@Home is in part powered by the <a href='http://boinc.berkeley.edu/'>Berkeley Open Infrastructure for Network Computing (BOINC)</a>.
                            </p>
                        </div>
                    </div>

             </div>
        </div> <!-- col-sm-6 -->

        </div> <!-- row -->
    </div> <!-- /container -->";



print_footer('Travis Desell and the DNA@Home Team', 'Travis Desell, Archana Dhasarathy, Sergei Nechaev');

echo "</body></html>";

?>
