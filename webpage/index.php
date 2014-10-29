<?php

$cwd[__FILE__] = __FILE__;
if (is_link($cwd[__FILE__])) $cwd[__FILE__] = readlink($cwd[__FILE__]);
$cwd[__FILE__] = dirname($cwd[__FILE__]);

require_once($cwd[__FILE__] . "/../../citizen_science_grid/header.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/navbar.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/news.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/footer.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/my_query.php");
require_once($cwd[__FILE__] . "/../../citizen_science_grid/csg_uotd.php");

print_header("DNA@Home", "", "dna");
print_navbar("Projects: DNA@Home", "DNA@Home", "..");

echo "
    <div class='container'>
        <div class='row'>
            <div class='col-sm-12'>
";


$carousel_info['items'][] = array(
    'active' => 'true',
    'image' => './images/dna_image_2.png',
    'text' => "DNA@Home is a joint effort between the <a href='http://cs.und.edu'>Computer Science Department</a> and the <a href='http://www.med.und.edu/basic-sciences/index.cfm'>Basic Sciences Department</a> of the <a href='http://und.edu'>University of North Dakota</a> and has been developed with support from <a href='http://rpi.edu'>Rensselaer Polytechnic Institute</a>. The goal of DNA@Home is to discover what regulates the genes in DNA. Ever notice that skin cells are different from a muscle cells, which are different from a bone cells, even though all these cells have every gene in your genome? That's because not all genes are \"on\" all the time. Depending on the cell type and what the cell is trying to do at any given moment, only a subset of the genes are used, and the remainder are shut off. DNA@home uses statistical algorithms to unlock the key to this differential regulation, using your volunteered computers.");

$carousel_info['items'][] = array(
    'image' => './images/dna-blue.png',
    'text' => 'The primary means by which genes are regulated is at the stage of "transcription" where a molecule called a polymerase reads along the DNA from the start of the gene to the end of the gene creating an RNA messenger. Other molecules, called transcription factors, bind to the DNA near the beginning of the gene and can help to recruit the polymerase or they can get in the way of, or inhibit, the polymerase. It is the presence or absence of the binding of these transcription factors that determine whether a gene is "on" or "off" but, for the most part, scientists do not know which transcription factors are responsible for regulating which genes.');

$carousel_info['items'][] = array(
    'image' => './images/protein_binding.png',
    'text' => "Transcription factors have \"fingers\" that prefer a certain short, sloppy pattern in the nucleotides \"letters\" of a DNA sequence, but in many cases we don't know what these patterns are. Our software looks for short sequences of nucleotides that appear more-or-less the same near multiple gene beginnings and which also appear more-or-less the same in the corresponding locations in the genomes of related species. As DNA sequences are huge, ranging from millions to billions of nucleotides, and these sequences are short and only approximately conserved from one site to the next, this is a real needle-in-the-haystack problem and requires lots of computational power. We hope that your computers can help.");

$carousel_info['items'][] = array(
    'image' => './images/e-cadherin_before_snail_expression.png',
    'caption' =>'This image shows what happens to the E-cadherin protein (stained in red) before Snail expression.',
    'text' => "<h4>Snail and <a href='http://en.wikipedia.org/wiki/Epithelial–mesenchymal_transition'>Epithelial–Mesenchymal Transition</a> (part 1):</h4><p>DNA@Home is currently investigating the <a href='http://en.wikipedia.org/wiki/SNAI1'>Snail</a> and <a href='http://en.wikipedia.org/wiki/SNAI2'>Slug'</a> transcription factors, which hav been shown to play important roles in both development and in disease. The Snail and Slug proteins are highly similar, with the repressive region in the first half of the proteins are 89% similar, while the end regions (the DNA-binding regions) are 84% identical, with most differences occuring the middle of the proteins. Understanding how these proteins know which specific DNA sequences to bind to is crucial in developing therapeutic targets for diseases like cancer.");

$carousel_info['items'][] = array(
    'image' => './images/e-cadherin_after_snail_expression.png',
    'caption' => 'This image shows what happens to the E-cadherin protein (stained in red) after Snail expression.',
    'text' => "<h4>Snail and <a href='http://en.wikipedia.org/wiki/Epithelial–mesenchymal_transition'>Epithelial–Mesenchymal Transition</a> (part 2):</h4><p> Expression of Snail and Slug causes loss of cell-adhesion proteins making the cells less \"sticky\" and more \"mobile\", causing metastasis. The cells are usually in an “epithelial” state, which means that the cells are sticking together. However, after Snail is expressed, you can see a nicely defined pattern of E-cadherin, which shows a cobblestone or chicken-wire pattern (red). Upon Snail induction, E-cadherin expression is lost, so the cells lose their “stickiness”, and start to move. The nucleus of each cell is stained in blue.</p>");

$projects_template = file_get_contents($cwd[__FILE__] . "/../../citizen_science_grid/templates/carousel.html");

$m = new Mustache_Engine;
echo $m->render($projects_template, $carousel_info);


echo "
            <div class='btn-group btn-group-justified' style='margin-top:20px;'>
                <a class='btn btn-primary' role='button' href='../csg/instructions.php'><h4>Volunteer Your Computer</h4></a>
            </div>

            </div> <!-- col-sm-12 -->
        </div> <!-- row -->

        <div class='row'>
            <div class='col-sm-6'>";

show_uotd(2, 10, "style='margin-top:20px;'");
csg_show_news();

echo "
            </div> <!-- col-sm-6 -->

            <div class='col-sm-6'>";

include $cwd[__FILE__] . "/templates/dna_info.html";

echo "
            </div> <!-- col-sm-6 -->
        </div> <!-- row -->
    </div> <!-- /container -->";


print_footer('Travis Desell and the DNA@Home Team', 'Travis Desell, Archana Dhasarathy, Sergei Nechaev');

echo "</body></html>";

?>
