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

// sample_work_generator: example BOINC work generator.
//
// --app name               app name (default example_app)
// --in_template_file       input template file (default example_app_in)
// --out_template_file      output template file (default example_app_out)
// -d N                     log verbosity level (0..4)
// --help                   show usage
// --version                show version
//
// - Runs as a daemon, and creates an unbounded supply of work.
//   It attempts to maintain a "cushion" of 100 unsent job instances
//   for the given app.
//   (your app may not work this way; e.g. you might create work in batches)
// - Creates a new input file for each job;
//   the file (and the workunit names) contain a timestamp
//   and sequence number, so they're unique.
//
// This is an example - customize for your needs

#include <sys/param.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>

#include <limits>

#include <fstream>
using std::ifstream;
using std::getline;

#include <sstream>
using std::ostringstream;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "backend_lib.h"
#include "boinc_db.h"
#include "error_numbers.h"
#include "filesys.h"
#include "parse.h"
#include "str_replace.h"
#include "str_util.h"
#include "svn_version.h"
#include "util.h"

#include "sched_config.h"
#include "sched_util.h"
#include "sched_msgs.h"

#include "mysql.h"

//from undvc_common
#include "file_io.hxx"

#define CUSHION 100
    // maintain at least this many unsent results
#define REPLICATION_FACTOR  2
    // number of instances of each job

const char* app_name = "gibbs";
const char* initial_workunit_xml_filename = "mtuberculosis_wu_initial_1.xml";
const char* checkpoint_workunit_xml_filename = "mtuberculosis_wu_checkpoint_1.xml";
const char* result_xml_filename = "mtuberculosis_result_1.xml";

char* in_template;
DB_APP app;
int start_time;
int seqno;

/**
 *  This wrapper makes for much more informative error messages when doing MYSQL queries
 */
#define mysql_query_check(conn, query) __mysql_check (conn, query, __FILE__, __LINE__)

void __mysql_check(MYSQL *conn, string query, const char *file, const int line) {
    mysql_query(conn, query.c_str());

    if (mysql_errno(conn) != 0) {
        ostringstream ex_msg;
        ex_msg << "ERROR in MySQL query: '" << query.c_str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << file << ":" << line;
        cerr << ex_msg.str() << endl;
        exit(1);
    }
}

int get_number_nucleotides(string sequence_filename) {
    ifstream sequence_file(sequence_filename.c_str());

    int number_nucleotides = 0;
    int line_number = 0;
    string line;
    while (getline(sequence_file, line)) {
        if (line.size() == 0) continue;
        if (line[0] == '>') continue;
        if (line.find_first_not_of(' ') == string::npos) continue;
        if (line.find_first_not_of("ACGTacgt") != string::npos) {
            cerr << "ERROR: malformed sequence file: '" << sequence_filename << "'" << endl;
            cerr << "Nucleotide that was not A, C, G, T, a c g or t." << endl;
            cerr << "Problem on line: " << line_number << endl;
            cerr << "line: '" << line << "'" << endl;
            continue;
        }

        number_nucleotides += line.size();
        line_number++;
    }
    cerr << "number_nucleotides: " << number_nucleotides << endl;

    return number_nucleotides;
}

// create one new job
int create_workunit(MYSQL *conn, const string &sampler_name, const string &sequences_filename,
                    int sampler_id, int wu_num, string command_line, double rsc_fpops_est, double rsc_fpops_bound,
                    double rsc_memory_bound, double rsc_disk_bound, double delay_bound) {

    DB_WORKUNIT wu;
    char name[256], path[MAXPATHLEN];
    const char* infiles[1];
    int retval;

    // make a unique name (for the job and its input file)
    ostringstream wu_name;
    wu_name << app_name << "_" << sampler_name << "_" << wu_num << "_0";

    sprintf(name, "%s", wu_name.str().c_str());   //should probably make this more descriptive

    // Fill in the job parameters
    wu.clear();
    wu.appid = app.id;
    safe_strcpy(wu.name, name);
    wu.rsc_fpops_est = rsc_fpops_est;
    wu.rsc_fpops_bound = rsc_fpops_bound;
    wu.rsc_memory_bound = rsc_memory_bound;
    wu.rsc_disk_bound = rsc_disk_bound; //50MB
    wu.delay_bound = delay_bound;
    wu.min_quorum = REPLICATION_FACTOR;
    wu.target_nresults = REPLICATION_FACTOR;
    wu.max_error_results = REPLICATION_FACTOR*4;
    wu.max_total_results = REPLICATION_FACTOR*8;
    wu.max_success_results = REPLICATION_FACTOR*4;
    wu.batch = sampler_id;

    //for initial workunits, the input file is only the sequencesFilename
    //for continuation workunits, there is the sequencesFilename and the final checkpoint of the previous workunit

    string only_filename = sequences_filename.substr(sequences_filename.find_last_of("/\\") + 1);
    cerr << "infile: '" << only_filename << "'" << endl;
    infiles[0] = only_filename.c_str();

    //need to make a gibbs_walk in the database and initialize these
    int current_steps = 0;  //current_steps is initially 0 as this is a new walk
    int seed = (int)(drand() * std::numeric_limits<int>::max());

    //random seed needs to be appended to the command line
    ostringstream extra_options;
    extra_options << " --seed " << seed;
    command_line.append(extra_options.str());

    //insert this gibbs_walk into the database
    ostringstream query;
    query   << "INSERT INTO gibbs_walk SET "
            << "seed = " << seed
            << ", sampler_id = " << sampler_id
            << ", current_steps = " << current_steps
            << ", current_sites = ''";
    mysql_query_check(conn, query.str().c_str());

    int walk_id = mysql_insert_id(conn);   //need to get this after inserting the walk into the database

    ostringstream additional_xml;
    additional_xml << "<seed>" << seed << "</seed>"
                   << "<sampler_id>" << sampler_id << "</sampler_id>"
                   << "<walk_id>" << walk_id << "</walk_id>"
                   << "<current_steps>" << current_steps << "</current_steps>";

    cerr << "additional_xml: '" << additional_xml.str() << "'" << endl;

    // Register the job with BOINC
    sprintf(path, "templates/%s", result_xml_filename);
    retval = create_work(
        wu,
        in_template,                    //input template file
        path,                           //result template file path
        config.project_path(path),
        infiles,                        //input files
        1,                              //number of input files
        config,
        command_line.c_str(),
        additional_xml.str().c_str()
    );

    return retval;
}

void main_loop(const vector<string> &arguments, MYSQL *conn) {
    string sampler_name = arguments[1];
    string sequences_filename = arguments[2];
    int number_walks = atoi(arguments[3].c_str());

    int workunit_steps   = 10000;
    int burn_in          = 0;                   //burn in for the full walk
    int samples          = 100000;             //samples to be taken for the full walk

    //check to see if the server is stopped
    check_stop_daemons();

    int numberNucleotides = get_number_nucleotides(sequences_filename);

    int numberMotifs = 4;
    int modelWidth = 6;
    ostringstream motif_string;
    motif_string << "forward," << modelWidth << " reverse," << modelWidth << " forward," << modelWidth << " reverse," << modelWidth << " ";
    int maxSites = 4;

    //Make sure the sequences filename is in the download directory
    copy_file_to_download_dir(sequences_filename);

    double rsc_fpops_est = ((double)(workunit_steps)) * (double)numberNucleotides * numberMotifs * (1.0 + (2.0 * modelWidth) + maxSites);

    cerr << ": " << workunit_steps << endl;
    cerr << ": " << (workunit_steps * numberNucleotides) << endl;
    cerr << ": " << (workunit_steps * numberNucleotides * numberMotifs) << endl;
    cerr << ": " << (workunit_steps * numberNucleotides * numberMotifs * (1.0 + (2.0 * modelWidth) + maxSites)) << endl;

    rsc_fpops_est *= 10.0;
    double rsc_fpops_bound = rsc_fpops_est * 100.0;
    double rsc_memory_bound = 5e7;
    double rsc_disk_bound = 50 * 1024 * 1024; //50MB
    double delay_bound = 86400;

    cerr << "NEW RSC_FPOPS_EST = " << rsc_fpops_est << endl;

    int workunit_burn_in = workunit_steps;
    int workunit_samples = 0;
    if (workunit_burn_in > burn_in) {
        workunit_burn_in = burn_in;
        workunit_samples = (workunit_steps - workunit_burn_in);
    }

    ostringstream command_line;
    command_line << " --max_sites " << maxSites
                 << " --blocks  0.1 0.225 0.225 0.225 0.225"
                 << " --motifs " + motif_string.str()
                 << " --enable_shifting 2 5"
                 << " --print_best_sites 0.1"
                 << " --print_current_sites"
//                 << " --print_accumulated_samples"
                 << " --checkpoint_frequency 1000"
                 << " --sequence_file sequences.txt"
                 << " --burn_in_period " << workunit_burn_in
                 << " --sample_period " << workunit_samples;


    //insert the new gibbs sampler into the database
    ostringstream query;
    query   << "INSERT INTO gibbs_sampler SET name = '" << sampler_name << "'"
            << ", initial_workunit_xml_filename = '" << initial_workunit_xml_filename << "'" 
            << ", checkpoint_workunit_xml_filename = '" << checkpoint_workunit_xml_filename << "'"
            << ", result_xml_filename = '" << result_xml_filename << "'"
            << ", sequences_filename = '" << sequences_filename << "'"
            << ", number_walks = " << number_walks
            << ", burn_in = " << burn_in
            << ", samples = " << samples
            << ", workunit_steps = " << workunit_steps
            << ", command_line_options = '" << command_line.str().c_str() << "'"
            << ", rsc_fpops_est = " << rsc_fpops_est
            << ", rsc_fpops_bound = " << rsc_fpops_bound
            << ", rsc_disk_bound = " << rsc_disk_bound
            << ", rsc_memory_bound = " << rsc_memory_bound
            << ", delay_bound = " << delay_bound;

    mysql_query_check(conn, query.str().c_str());
    int sampler_id = mysql_insert_id(conn); //get this from inserting the gibbs_walk into the database

    for (int i = 0; i < number_walks; i++) {
        create_workunit(conn, sampler_name, sequences_filename, sampler_id, i, command_line.str(), rsc_fpops_est, rsc_fpops_bound, rsc_memory_bound, rsc_disk_bound, delay_bound);
    }
}

void usage(char *name) {
    fprintf(stderr, "This is an example BOINC work generator.\n"
        "This work generator has the following properties\n"
        "(you may need to change some or all of these):\n"
        "  It attempts to maintain a \"cushion\" of 100 unsent job instances.\n"
        "  (your app may not work this way; e.g. you might create work in batches)\n"
        "- Creates work for the application \"example_app\".\n"
        "- Creates a new input file for each job;\n"
        "  the file (and the workunit names) contain a timestamp\n"
        "  and sequence number, so that they're unique.\n\n"
        "Usage: %s [OPTION]...\n\n"
        "Options:\n"
        "  [ --app X                Application name (default: example_app)\n"
        "  [ --in_template_file     Input template (default: example_app_in)\n"
        "  [ --out_template_file    Output template (default: example_app_out)\n"
        "  [ -d X ]                 Sets debug level to X.\n"
        "  [ -h | --help ]          Shows this help text.\n"
        "  [ -v | --version ]       Shows version information.\n",
        name
    );
}

int main(int argc, char** argv) {
    int i, retval;
    char buf[256];

    for (i=1; i<argc; i++) {
        if (is_arg(argv[i], "d")) {
            if (!argv[++i]) {
                log_messages.printf(MSG_CRITICAL, "%s requires an argument\n\n", argv[--i]);
                usage(argv[0]);
                exit(1);
            }
            int dl = atoi(argv[i]);
            log_messages.set_debug_level(dl);
            if (dl == 4) g_print_queries = true;
        } else if (is_arg(argv[i], "h") || is_arg(argv[i], "help")) {
            usage(argv[0]);
            exit(0);
        } else if (is_arg(argv[i], "v") || is_arg(argv[i], "version")) {
            printf("%s\n", SVN_VERSION);
            exit(0);
//        } else {
//            log_messages.printf(MSG_CRITICAL, "unknown command line argument: %s\n\n", argv[i]);
//            usage(argv[0]);
//            exit(1);
        }
    }

    //Aaron Comment: if at any time the retval value is greater than 0, then the program
    //has failed in some manner, and the program then exits.

    //Aaron Comment: processing project's config file.
    retval = config.parse_file();
    if (retval) {
        log_messages.printf(MSG_CRITICAL,
            "Can't parse config.xml: %s\n", boincerror(retval)
        );
        exit(1);
    }

    //Aaron Comment: opening connection to database.
    retval = boinc_db.open(config.db_name, config.db_host, config.db_user, config.db_passwd);
    if (retval) {
        log_messages.printf(MSG_CRITICAL, "can't open db\n");
        exit(1);
    }

    //Aaron Comment: looks for applicaiton to be run. If not found, program exits.
    sprintf(buf, "where name='%s'", app_name);
    if (app.lookup(buf)) {
        log_messages.printf(MSG_CRITICAL, "can't find app %s\n", app_name);
        exit(1);
    }

    //Aaron Comment: looks for work templates, if cannot find, or are corrupted,
    //the program exits.
    sprintf(buf, "templates/%s", initial_workunit_xml_filename);
    if (read_file_malloc(config.project_path(buf), in_template)) {
        log_messages.printf(MSG_CRITICAL, "can't read input template %s\n", buf);
        exit(1);
    }

    //Aaron Comment: if work generator passes all startup tests, the main work gneration
    //loop is called.
    start_time = time(0);
    seqno = 0;

    log_messages.printf(MSG_NORMAL, "Starting\n");

    main_loop(vector<string>(argv, argv + argc), boinc_db.mysql);
}
