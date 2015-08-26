<?php

foreach (glob("*.out") as $filename) {
    $command = "sed -i -e 1,6d $filename";
    echo $command . "\n";
    exec($command);
}


?>
