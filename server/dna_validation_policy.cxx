// This file is part of BOINC.
// http://boinc.berkeley.edu
// Copyright (C) 2008 University of California
//
// BOINC is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// BOINC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with BOINC.  If not, see <http://www.gnu.org/licenses/>.

// A sample validator that requires a majority of results to be
// bitwise identical.
// This is useful only if either
// 1) your application does no floating-point math, or
// 2) you use homogeneous redundancy

#include "config.h"
#include "util.h"
#include "sched_util.h"
#include "sched_msgs.h"
#include "validate_util.h"
#include "md5_file.h"
#include "error_numbers.h"
#include "stdint.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>


#include "undvc_common/parse_xml.hxx"
#include "undvc_common/file_io.hxx"

using std::string;
using std::vector;
using std::ifstream;

using std::string;
using std::endl;

int init_result(RESULT& result, void*& data) {
    char *probabilities;

    try {
        string prob_str = parse_xml<string>(result.stderr_out, "slice_probabilities");

        replace( prob_str.begin(), prob_str.end(), '\n', ' ' );

        char chars[] = "\n\r";
        for (unsigned int i = 0; i < strlen(chars); ++i) {
            // you need include <algorithm> to use general algorithms like std::remove()
            prob_str.erase (std::remove(prob_str.begin(), prob_str.end(), chars[i]), prob_str.end());
        }

        probabilities = (char*)malloc(sizeof(char) * (strlen(prob_str.c_str()) + 1));
        strcpy(probabilities, prob_str.c_str());
        probabilities[strlen(prob_str.c_str())] = '\0';

        cout << "probabilities: " << probabilities << endl;
    } catch (string error_message) {
        log_messages.printf(MSG_CRITICAL, "wildlife_validation_policy get_data_from_result([RESULT#%d %s]) failed with error: %s\n", result.id, result.name, error_message.c_str());
        log_messages.printf(MSG_CRITICAL, "XML:\n%s\n", result.stderr_out);
        result.outcome = RESULT_OUTCOME_VALIDATE_ERROR;
        result.validate_state = VALIDATE_STATE_INVALID;

//        exit(1);
        return ERR_XML_PARSE;
    }

    data = (void*)probabilities;

    return 0;
}

int compare_results(
    RESULT & r1, void* data1,
    RESULT const& r2, void* data2,
    bool& match
) {
    char *probabilities1 = (char*)data1;
    char *probabilities2 = (char*)data2;

    vector<double> p1;
    istringstream iss1(probabilities1);
    copy(istream_iterator<double>(iss1), istream_iterator<double>(), back_inserter<vector<double> >(p1));
    
    vector<double> p2;
    istringstream iss2(probabilities2);
    copy(istream_iterator<double>(iss2), istream_iterator<double>(), back_inserter<vector<double> >(p2));

    if (p1.size() != p2.size()) {
        match = false;
        log_messages.printf(MSG_CRITICAL, "ERROR, number of probabilities is different. %d vs %d\n", (int)p1.size(), (int)p2.size());

        /*
        log_messages.printf(MSG_CRITICAL, "p1 string: '%s'\n", probabilities1);
        log_messages.printf(MSG_CRITICAL, "p2 string: '%s'\n", probabilities2);

        log_messages.printf(MSG_CRITICAL, "probabilities1:\n");
        for (uint32_t i = 0; i < p1.size(); i++) {
            log_messages.printf(MSG_CRITICAL, "\t%lf\n", p1[i]);
        }

        log_messages.printf(MSG_CRITICAL, "probabilities2:\n");
        for (uint32_t i = 0; i < p2.size(); i++) {
            log_messages.printf(MSG_CRITICAL, "\t%lf\n", p2[i]);
        }
        */

        match = false;

        return 0;
    }

    double threshold = 0.026;

    for (uint32_t i = 0; i < p1.size(); i++) {
        if (fabs(p1[i] - p2[i]) > threshold) {
            match = false;

            /*
            log_messages.printf(MSG_CRITICAL, "probabilities1:\n");
            for (uint32_t j = 0; j < p1.size(); j++) {
                log_messages.printf(MSG_CRITICAL, "\t%lf\n", p1[j]);
            }

            log_messages.printf(MSG_CRITICAL, "probabilities2:\n");
            for (uint32_t j = 0; j < p2.size(); j++) {
                log_messages.printf(MSG_CRITICAL, "\t%lf\n", p2[j]);
            }
            */

            log_messages.printf(MSG_CRITICAL, "ERROR, difference in probabilities (%lf) exceeded threshold (%lf):\n", fabs(p1[i]-p2[i]), threshold);
            log_messages.printf(MSG_CRITICAL, "probabilities1[%d]: %lf\n", i, p1[i]);
            log_messages.printf(MSG_CRITICAL, "probabilities2[%d]: %lf\n", i, p2[i]);
            exit(1);

            return 0;
        }
    }

    match = true;

    return 0;
}

int cleanup_result(RESULT const& /*result*/, void* data) {
    char* result = (char*)data;

    delete result;

    return 0;
}

const char *BOINC_RCSID_7ab2b7189c = "$Id: sample_bitwise_validator.cpp 21735 2010-06-12 22:08:15Z davea $";
