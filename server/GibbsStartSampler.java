package gibbs_assimilator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;

import java.util.Arrays;

import jboinc.util.DirectoryTree;
import jboinc.util.Parameters;

import jfgdo.GibbsSampler;

public class GibbsStartSampler {

	private static String initial_workunit_xml_filename = "templates/mtuberculosis_wu_initial_1.xml";
	private static String checkpoint_workunit_xml_filename = "templates/mtuberculosis_wu_checkpoint_1.xml";
	private static String burn_in_result_xml_filename = "templates/mtuberculosis_burn_in_result_1.xml";
	private static String sampling_result_xml_filename = "templates/mtuberculosis_sampling_result_1.xml";

	public static void main(String[] arguments) {
        if (arguments.length != 4) {
            System.err.println("Error running GibbsStartSampler, usage: ");
            System.err.println("\tjava GibbsStartSampler <sampler name> <sequences filename> <workunitBurnIn> <samplesBurnIn>");
            System.exit(0);
        }

		DirectoryTree.setBaseDirectory("/projects/dna/");
		DirectoryTree.setBaseURL("http://volunteer.cs.und.edu/dna/");
		DirectoryTree.setResultsDirectory("/projects/dna/results/");

		String samplerName = arguments[0];
		String sequencesFilename = arguments[1];
		File sequencesFile = new File( sequencesFilename );

        int workunitBurnIn = Integer.parseInt(arguments[2]);
        int workunitSamples = Integer.parseInt(arguments[3]);

        try {
                System.out.println("copying sequencesFilename: " + sequencesFile + ", getName(): " + sequencesFile.getName());
                DirectoryTree.copyToDownloadDirectory( sequencesFile );

        } catch (Exception e) {
                System.err.println("Could not copy files to download directory, exception thrown: " + e);
                e.printStackTrace();
                return;
        }
        sequencesFilename = sequencesFile.getName();

        int numberNucleotides = getNumberNucleotides(sequencesFile);

        int numberMotifs = 2;
        int modelWidth = 16;
        String motif_string = "forward," + modelWidth + " reverse," + modelWidth + " ";
        int maxSites = 3;

		int numberWalks = 10000;
		int burnIn = 100000;
		int samples = 100000;
//		int workunitBurnIn = 100000;
//		int workunitSamples = 100000;

        double rsc_fpops_est = ((double)(workunitBurnIn + workunitSamples)) * (double)numberNucleotides * numberMotifs * (1.0 + (2.0 * modelWidth) + maxSites);

        System.err.println(": " + (workunitBurnIn + workunitSamples));
        System.err.println(": " + ((workunitBurnIn + workunitSamples) * numberNucleotides) );
        System.err.println(": " + ((workunitBurnIn + workunitSamples) * numberNucleotides * numberMotifs) );
        System.err.println(": " + ((workunitBurnIn + workunitSamples) * numberNucleotides * numberMotifs * (1.0 + (2.0 * modelWidth) + maxSites)) );

        rsc_fpops_est *= 10.0; 
        double rsc_fpops_bound = rsc_fpops_est * 100.0;

        System.err.println("NEW RSC_FPOPS_EST = " + rsc_fpops_est);

        System.err.println("OLD RSC_FPOPS_EST = " + 4 * 1024.0 * 1024.0 * 1024.0 * 1024.0);
//        System.exit(0);

		String commandLineOptions =
                                " --max_sites " + maxSites + " " +
                                " --blocks 0.1 0.3 0.3 0.3 " +
                                " --motifs " + motif_string +
                                " --enable_shifting 2 5 " +
                                " --print_best_sites 0.1 " +
                                " --print_current_sites "+
                                " --print_accumulated_samples ";

		new GibbsSampler(samplerName, initial_workunit_xml_filename, checkpoint_workunit_xml_filename, burn_in_result_xml_filename, sampling_result_xml_filename, sequencesFilename, numberWalks, burnIn, samples, workunitBurnIn, workunitSamples, commandLineOptions, rsc_fpops_bound, rsc_fpops_est);
	}

    public static int getNumberNucleotides(File sequencesFile) {
        int numberNucleotides = 0;

        try {
            BufferedReader in = new BufferedReader( new FileReader( sequencesFile ) );

            String line = in.readLine();
            while (line != null) {
                if ( !(line.length() == 0 || line.charAt(0) == '>') ) {
                    numberNucleotides += line.length();

                    if ( !Character.isLetter(line.charAt(line.length() - 1)) ) {
                        numberNucleotides--;
                    }
                }
                line = in.readLine();
//                System.err.println(line);
            }
            in.close();

        } catch (Exception e) {
            System.err.println("Error reading sequences file: " + e);
            e.printStackTrace();
            System.exit(0);
        }

        System.out.println("number of nucleotides: " + numberNucleotides);

        return numberNucleotides;
    }
}
