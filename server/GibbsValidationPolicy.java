package gibbs_assimilator;

import jboinc.daemons.ValidationPolicy;

import jboinc.database.Host;
import jboinc.database.Result;
import jboinc.database.Workunit;

import jfgdo.GibbsSampler;

import jboinc.util.DirectoryTree;
import jboinc.util.XMLTemplate;
import jboinc.util.XMLParseException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedList;
import java.text.NumberFormat;

public class GibbsValidationPolicy implements ValidationPolicy {
	public final static double MAXIMUM_ERROR = 10e-10;

	private NumberFormat nf = NumberFormat.getInstance();
	public GibbsValidationPolicy() {
		nf.setMinimumFractionDigits(20);
	}

	public boolean checkPair(Result canonicalResult, Result unvalidatedResult) {
		String unvalidatedResultString = null;

		String resultName = unvalidatedResult.getName();
		String searchName = resultName.substring(0, resultName.lastIndexOf('_', resultName.lastIndexOf('_') - 1) - 1);

		try {
			unvalidatedResultString = XMLTemplate.processString(unvalidatedResult.getStderrOut(), "current_sites");

		} catch (XMLParseException e) {
			System.err.println("Could not get unvalidated search result, xml parse exception.");
			System.err.println("Exception: " + e);
			e.printStackTrace();

			unvalidatedResult.setOutcome(Result.RESULT_OUTCOME_VALIDATE_ERROR);
			return false;
		}

		String canonicalResultString = null;
		try {
			canonicalResultString = XMLTemplate.processString(canonicalResult.getStderrOut(), "current_sites");

		} catch (XMLParseException e) {
			System.err.println("Could not get unvalidated search result, xml parse exception.");
			System.err.println("Exception: " + e);
			e.printStackTrace();

			unvalidatedResult.setOutcome(Result.RESULT_OUTCOME_VALIDATE_ERROR);
			System.err.println("ERROR PARSING CANONICAL RESULT XML!\n");
			System.err.println("canonicalResult: " + canonicalResult);
			System.exit(0);
		}


		if (unvalidatedResultString.equals(canonicalResultString)) {
			System.out.println("result: " + resultName + " is valid.");

			unvalidatedResult.setValidateState(Result.VALIDATE_STATE_VALID);
		} else {
			unvalidatedResult.setValidateState(Result.VALIDATE_STATE_INVALID);
		}

		return false;
	}

	public boolean checkSet(Workunit workunit, LinkedList<Result> unvalidatedResults) {
		int i, j;

		String workunit_template = workunit.getWorkunitTemplate();

		String[] resultSites = new String[unvalidatedResults.size()];

		i = 0;
		for (Result unvalidatedResult : unvalidatedResults) {
			try {
				resultSites[i] = XMLTemplate.processString(unvalidatedResult.getStderrOut(), "current_sites");

			} catch (XMLParseException e) {
				System.err.println("Could not get unvalidated search result, xml parse exception.");
				System.err.println("Exception: " + e);
				e.printStackTrace();

				unvalidatedResult.setOutcome(Result.RESULT_OUTCOME_VALIDATE_ERROR);
				resultSites[i] = "";
			}

			i++;
		}


		int quorum = workunit.getMinQuorum();
		if (quorum < 2) {
			System.out.println("ERROR: workunit quorum < 2!");
			System.out.println(workunit);
			System.exit(0);
		}

		for (i = 0; i < resultSites.length; i++) {
			System.out.println("resultSites[" + i + "].length() before replaceAll: " + resultSites[i].length());
			resultSites[i] = resultSites[i].replaceAll("\r", "");
			System.out.println("resultSites[" + i + "].length() after  replaceAll: " + resultSites[i].length());
		}

		int similarResults = 0;
		Result canonicalResult = null;
		String canonicalResultString = null;
		for (i = 0; i < resultSites.length; i++) {
			if (resultSites[i].equals("")) continue;

			similarResults = 1;
			for (j = i + 1; j < resultSites.length; j++) {
				if (resultSites[i].equals(resultSites[j])) similarResults++;
			}

			if (similarResults >= quorum) {
				canonicalResultString = resultSites[i];
				canonicalResult = unvalidatedResults.get(i);
				workunit.setCanonicalResultId(canonicalResult.getId());

				GibbsSampler.insertWorkunit(workunit, canonicalResultString);

				System.out.println("set canonical result: " + canonicalResult);
				break;
			}
		}	

		if (canonicalResult == null) return false;

		for (i = 0; i < resultSites.length; i++) {
			if (resultSites[i].equals(canonicalResultString)) {
				unvalidatedResults.get(i).setValidateState(Result.VALIDATE_STATE_VALID);
				System.out.println("set resultSites[" + i + "] to valid.");
			} else {
				unvalidatedResults.get(i).setValidateState(Result.VALIDATE_STATE_INVALID);
				System.out.println("set resultSites[" + i + "] to invalid.");
			}
		}

		return false;
	}

}
