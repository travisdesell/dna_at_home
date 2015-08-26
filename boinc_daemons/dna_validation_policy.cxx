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

char* read_file(OUTPUT_FILE_INFO &fi) {
    ifstream samples_file(fi.path.c_str());

    string file_contents, line;
    while (getline(samples_file, line)) {
        file_contents.append(line);
    }

    file_contents.erase(std::remove(file_contents.begin(), file_contents.end(), '\r'), file_contents.end());
    file_contents.erase(std::remove(file_contents.begin(), file_contents.end(), '\n'), file_contents.end());

    char* contents = new char[file_contents.size() + 1];

    //sprintf appends the null character
    strcpy(contents, file_contents.c_str());
    contents[file_contents.size()] = '\0';

    return contents;
}

int init_result(RESULT& result, void*& data) {
    ResultData *rd = new ResultData[1];

    //need to see if current samples are in a file
    vector<OUTPUT_FILE_INFO> files;

    int retval = get_output_file_infos(result, files);
    if (retval) {
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: can't get output filenames\n", result.id, result.name);  
        return retval;
    }   

    if (files.size() > 2) {
        //some error, there should be only one or two files.
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: too many output filenames\n", result.id, result.name);  
        result.outcome = RESULT_OUTCOME_VALIDATE_ERROR;
        result.validate_state = VALIDATE_STATE_INVALID;
        exit(1);
        return 1;

    } else if (files.size() == 2) {
        rd->current_sites = read_file(files[0]);
//        cout << "sites: " << rd->current_sites << endl << endl;
        rd->current_samples = read_file(files[1]);
//        cout << "samples: " << rd->current_samples << endl << endl;

    } else if (files.size() == 1) {
        //theres one file, it will be the samples
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: need to read samples from output file\n", result.id, result.name);  

        rd->current_sites = read_file(files[0]);
//        cout << "sites: " << rd->current_sites << endl << endl;
//        cout << "samples: NULL" << endl << endl;

        rd->current_samples = NULL;

        exit(1);

    } else {
        //no files, samples aren't being recorded yet
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] check_set: not enough output filenames\n", result.id, result.name);  
        result.outcome = RESULT_OUTCOME_VALIDATE_ERROR;
        result.validate_state = VALIDATE_STATE_INVALID;
        exit(1);
        return 1;
    }

    data = (void*)rd;

    return 0;
}

void split_samples(const string &s, vector< vector<int> > &values) {
    values.clear();

    stringstream line_stream(s);
    string line;
    int current = 0;
    //cout << "parsing samples!" << endl;
    while (getline(line_stream, line, '\n')) {
        stringstream value_stream(line);
        string value;

        values.push_back(vector<int>());

        while (getline(value_stream, value, ',')) {
            values[current].push_back( atoi(value.c_str()) );

            //cout << " " << values[current].back();
        }
        current++;
        //cout << endl;
    }
}

int compare_results(
    RESULT & r1, void* data1,
    RESULT const& r2, void* data2,
    bool& match
) {
    cerr << "COMPARING RESULTS!" << endl;

    ResultData *rd1 = (ResultData*)data1;
    ResultData *rd2 = (ResultData*)data2;

    if (strcmp(rd1->current_samples, rd2->current_samples) == 0) {
        if (strcmp(rd1->current_sites, rd2->current_sites) == 0) {
            log_messages.printf(MSG_CRITICAL, "sites and samples match!\n");
            match = true;
            return 0;

        } else {
            log_messages.printf(MSG_CRITICAL, "ERROR, current_sites are different.\n");
            //log_messages.printf(MSG_CRITICAL, "ERROR, current_sites are different.\n%s\nvs\n%s\n", rd1->current_sites, rd2->current_sites);

            match = false;
            //exit(1);
            return 0;
        }

    } else {
        log_messages.printf(MSG_CRITICAL, "ERROR, current_samples are different.\n");
        //log_messages.printf(MSG_CRITICAL, "ERROR, current_samples are different.\n%s\nvs\n%s\n", rd1->current_samples, rd2->current_samples);

        vector< vector<int> > samples1;
        vector< vector<int> > samples2;

        split_samples(rd1->current_samples, samples1);
        split_samples(rd2->current_samples, samples2);

        log_messages.printf(MSG_CRITICAL, "samples1 size: %zu, samples2 size: %zu\n", samples1.size(), samples2.size());

        if (samples1.size() != samples2.size()) {
            log_messages.printf(MSG_CRITICAL, "ERROR, samples1 and samples2 have a different number of sequences: %zu vs %zu.\n", samples1.size(), samples2.size());
            //exit(1);
            match = false;
            return 0;
        }

        uint64_t sum_difference = 0, total_sites = 0;
        uint64_t sum1 = 0, sum2 = 0;
        uint64_t max_difference = 0;
        uint64_t diff = 0;
        uint64_t diffs_over_100 = 0;
        for (uint32_t i = 0; i < samples1.size(); i++) {
            log_messages.printf(MSG_CRITICAL, "samples1[%d] size: %zu, samples2[%d] size: %zu\n", i, samples1[i].size(), i, samples2[i].size());

            if (samples1[i].size() != samples2[i].size()) {
                log_messages.printf(MSG_CRITICAL, "ERROR, samples1[%d] and samples2[%d] have a different number of nucleotides: %zu vs %zu.\n", i, i, samples1[i].size(), samples2[i].size());
                //exit(1);
            }

            for (uint32_t j = 0; j < samples1[i].size(); j++) {
                sum1 += samples1[i][j];
                sum2 += samples2[i][j];
                diff = (uint64_t)abs(samples1[i][j] - samples2[i][j]);
                sum_difference += diff;
                total_sites++;

                if (max_difference < diff) max_difference = diff;
                if (diff > 100) diffs_over_100++;
            }
        }
        log_messages.printf(MSG_CRITICAL, "sum of samples1: %ld, sum of samples2: %ld, max difference: %ld, diffs_over_100: %ld\n", sum1, sum2, max_difference, diffs_over_100);
        log_messages.printf(MSG_CRITICAL, "ERROR, sum difference of samples was: %ld, total sites sampled: %ld, average difference: %lf\n", sum_difference, total_sites, (double)sum_difference/(double)total_sites);

        match = false;
        //exit(1);
        return 0;
    }

    return 0;
}

int cleanup_result(RESULT const& /*result*/, void* data) {
    ResultData* rd = (ResultData*)data;

    if (rd->current_samples != NULL) delete[] rd->current_samples;
    if (rd->current_sites != NULL) delete[] rd->current_sites;
    delete rd;

    return 0;
}

const char *BOINC_RCSID_7ab2b7189c = "$Id: sample_bitwise_validator.cpp 21735 2010-06-12 22:08:15Z davea $";
