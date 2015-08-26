<?php

foreach (glob("*.out") as $filename) {
    $outfile = str_replace("out", "png", $filename);

    $command = "python plot_distances2.py $filename $outfile";

    echo $command . "\n";
    exec($command);
}

?>
