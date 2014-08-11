package gibbs_assimilator;

import jboinc.daemons.ValidatorNew;
import jboinc.util.DirectoryTree;

public class GibbsValidator {

	public static void main(String[] arguments) {
		DirectoryTree.setBaseDirectory("/export/www/dnahome/");
		DirectoryTree.setBaseURL("http://dnahome.cs.rpi.edu/dna/");
		DirectoryTree.setResultsDirectory("/export/www/dnahome/results/");

		ValidatorNew validator = new ValidatorNew(1, 0, Integer.parseInt(arguments[0]), new GibbsValidationPolicy());
		validator.start();
	}
}
