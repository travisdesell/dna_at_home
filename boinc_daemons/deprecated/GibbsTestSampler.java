package gibbs_assimilator;

import java.io.File;
import java.io.IOException;

import java.util.Arrays;

import jboinc.util.DirectoryTree;
import jboinc.util.Parameters;

import jfgdo.GibbsSampler;

public class GibbsTestSampler {

	private static String initial_workunit_xml_filename = "templates/mtuberculosis_wu_initial_1.xml";
	private static String checkpoint_workunit_xml_filename = "templates/mtuberculosis_wu_checkpoint_1.xml";
	private static String result_xml_filename = "templates/mtuberculosis_result_1.xml";

	public static void main(String[] arguments) {
		DirectoryTree.setBaseDirectory("/export/www/dnahome/");
		DirectoryTree.setBaseURL("http://dnahome.cs.rpi.edu/dna/");
		DirectoryTree.setResultsDirectory("/export/www/results/");

		new GibbsSampler(Integer.parseInt(arguments[0]));
	}
}
