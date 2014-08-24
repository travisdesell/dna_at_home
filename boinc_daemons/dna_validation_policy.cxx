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


//from undvc_common
#include "parse_xml.hxx"
#include "file_io.hxx"

using std::string;
using std::vector;
using std::ifstream;

using std::string;
using std::endl;

struct ResultData {
    char *current_sites;
    char *current_samples;
};

int init_result(RESULT& result, void*& data) {
    ResultData *rd = new ResultData[1];

    //need to see if current samples are in a file
    vector<OUTPUT_FILE_INFO> files;

    int retval = get_output_file_infos(result, files);
    if (retval) {
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: can't get output filenames\n", result.id, result.name);  
        return retval;
    }   

    if (files.size() > 1) {
        //some error, there should only be one file.
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: too many output filenames\n", result.id, result.name);  
        exit(1);
        return 1;

    } else if (files.size() == 1) {
        //theres one file, it will be the samples
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: need to read samples from output file\n", result.id, result.name);  

        OUTPUT_FILE_INFO& fi = files[0];
        ifstream samples_file(fi.path.c_str());

        string current_samples, line;
        while (getline(samples_file, line)) {
            current_samples.append(line);
        }

        current_samples.erase(std::remove(current_samples.begin(), current_samples.end(), '\r'), current_samples.end());
        current_samples.erase(std::remove(current_samples.begin(), current_samples.end(), '\n'), current_samples.end());

        rd->current_samples = new char[current_samples.size() + 1];

        //sprintf appends the null character
        strcpy(rd->current_samples, current_samples.c_str());

        rd->current_samples[current_samples.size()] = '\0';

    } else {
        //no files, samples aren't being recorded yet
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: samples aren't recorded yet\n", result.id, result.name);  

        rd->current_samples = (char*)malloc(sizeof(char) * 1);
        rd->current_samples[0] = '\0';
    }

    try {
        //need to get current_sites from XML
        //need to get current_samples from file (if exists) -- maybe get it from XML as well?

        //need to copy stderr_out into a string otherwise we get
        //segmentation faults for reusing that data
        string stderr_out(result.stderr_out);

        string current_sites = parse_xml<string>(stderr_out, "current_sites");

        //remove the \r from windows results
        //you need include <algorithm> to use general algorithms like std::remove()
        current_sites.erase(std::remove(current_sites.begin(), current_sites.end(), '\r'), current_sites.end());
        current_sites.erase(std::remove(current_sites.begin(), current_sites.end(), '\n'), current_sites.end());

        rd->current_sites = new char[current_sites.size() + 1];

        //sprintf appends the null character
        strcpy(rd->current_sites, current_sites.c_str());

        rd->current_sites[current_sites.size()] = '\0';

    } catch (string error_message) {
        log_messages.printf(MSG_CRITICAL, "dna_validation_policy get_data_from_result([RESULT#%d %s]) failed with error: %s\n", result.id, result.name, error_message.c_str());
        log_messages.printf(MSG_CRITICAL, "XML:\n%s\n", result.stderr_out);
        result.outcome = RESULT_OUTCOME_VALIDATE_ERROR;
        result.validate_state = VALIDATE_STATE_INVALID;

        rd->current_sites = new char[1];
        rd->current_sites[0] = '\0';
        data = (void*)rd;

        return ERR_XML_PARSE;
    }

    data = (void*)rd;

    return 0;
}

int compare_results(
    RESULT & r1, void* data1,
    RESULT const& r2, void* data2,
    bool& match
) {
    cerr << "COMPARING RESULTS!" << endl;

    ResultData *rd1 = (ResultData*)data1;
    ResultData *rd2 = (ResultData*)data2;

    if (strcmp(rd1->current_sites, rd2->current_sites) == 0) {
        if (strcmp(rd1->current_samples, rd2->current_samples) == 0) {
            log_messages.printf(MSG_CRITICAL, "sites and samples match!\n");
            match = true;
            return 0;

        } else {
            log_messages.printf(MSG_CRITICAL, "ERROR, current_samples are different.\n%s\nvs\n%s\n", rd1->current_samples, rd2->current_samples);
            match = false;
            return 0;
        }

    } else {
        log_messages.printf(MSG_CRITICAL, "ERROR, current_sites are different.\n%s\nvs\n%s\n", rd1->current_sites, rd2->current_sites);
        match = false;
        return 0;
    }

    return 0;
}

int cleanup_result(RESULT const& /*result*/, void* data) {
    ResultData* rd = (ResultData*)data;

    delete[] rd->current_samples;
    delete[] rd->current_sites;
    delete rd;

    return 0;
}

const char *BOINC_RCSID_7ab2b7189c = "$Id: sample_bitwise_validator.cpp 21735 2010-06-12 22:08:15Z davea $";
